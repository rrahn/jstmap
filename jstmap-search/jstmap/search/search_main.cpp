// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <functional>
#include <numeric>
#include <omp.h>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/range/to.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/all_matches.hpp>
#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/bam_writer.hpp>
#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/search_matches.hpp>
#include <jstmap/search/filter_queries.hpp>
#include <jstmap/search/match_aligner.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_main.hpp>
#include <jstmap/search/bucket_searcher.hpp>
#include <jstmap/search/bucket.hpp>
#include <jstmap/search/type_alias.hpp>
#include <jstmap/search/options.hpp>

namespace jstmap
{

int search_main(seqan3::argument_parser & search_parser)
{
    search_options options{};

    search_parser.add_positional_option(options.jst_input_file_path,
                                       "The path to the journaled sequence tree.",
                                       seqan3::input_file_validator{{"jst"}});
    search_parser.add_positional_option(options.query_input_file_path,
                                       "The path to the read file.",
                                       seqan3::input_file_validator{{"fa", "fasta"}});
    search_parser.add_positional_option(options.map_output_file_path,
                                       "The alignment map output file.",
                                       seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                     {"sam", "bam"}});

    search_parser.add_flag(options.is_quite,
                           'q',
                           "quite",
                           "Disables all logging.",
                           seqan3::option_spec::standard);
    search_parser.add_flag(options.is_verbose,
                           'v',
                           "verbose",
                           "Enables expansive debug logging.",
                           seqan3::option_spec::standard);

    search_parser.add_option(options.index_input_file_path,
                             'i',
                             "index",
                             "The prebuilt index to speedup the search.",
                             seqan3::option_spec::standard,
                             seqan3::input_file_validator{{"ibf"}});
    search_parser.add_option(options.error_rate,
                             'e',
                             "error-rate",
                             "The error rate allowed for mapping the reads.",
                             seqan3::option_spec::standard,
                             seqan3::arithmetic_range_validator{0.0, 1.0});
    search_parser.add_option(options.thread_count,
                             't',
                             "thread-count",
                             "The number of threads to use for the search.",
                             seqan3::option_spec::standard,
                             seqan3::arithmetic_range_validator{1u, std::thread::hardware_concurrency()});

    try
    {
        search_parser.parse();
        if (options.is_quite) {
            get_application_logger().set_verbosity(verbosity_level::quite);
        } else if (options.is_verbose) {
            get_application_logger().set_verbosity(verbosity_level::verbose);
        }

        log_debug("References file:", options.jst_input_file_path.string());
        log_debug("Query file:", options.query_input_file_path.string());
        log_debug("Output file:", options.map_output_file_path.string());
        log_debug("Index file:", options.index_input_file_path.string());
        log_debug("Error rate:", options.error_rate);
        log_debug("Thread count:", options.thread_count);
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        log_err(ex.what());
        return -1;
    }

    // Load the queries and the jst.
    auto global_start = std::chrono::high_resolution_clock::now();
    try
    {
        log_info("Start mapping");
        log_debug("Load reads");
        auto start = std::chrono::high_resolution_clock::now();
        std::vector query_records = load_queries(options.query_input_file_path);
        size_t query_idx{};
        auto queries = query_records
                     | std::views::transform([&] (sequence_record_t & record) {
                        return search_query{query_idx++, std::move(record)};
                     })
                     | seqan3::ranges::to<std::vector>();
        auto end = std::chrono::high_resolution_clock::now();
        log_debug("Read count", queries.size());
        log_debug("Loading time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        log_debug("Load reference database");
        start = std::chrono::high_resolution_clock::now();
        rcs_store_t rcs_store = load_jst(options.jst_input_file_path);
        // For debug purposes!

        // rcs_store_t rcs_store{reference_t{rcs_store_tmp.source().begin(), rcs_store_tmp.source().end()}, rcs_store_tmp.size()};
        // auto it_end = rcs_store_tmp.variants().end() - 1;
        // auto it = it_end - 1;
        // for (; it != it_end; ++it) {
        //     rcs_store.add(*it);
        // }

        end = std::chrono::high_resolution_clock::now();
        log_info("Loading time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        // Now we want to handle the bins as well.

        start = std::chrono::high_resolution_clock::now();
        size_t bin_size{std::numeric_limits<size_t>::max()};
        std::vector<search_queries_type> search_queries{};
        // What if the ibf is not present?
        if (options.index_input_file_path.empty())
        {
            log_debug("No prefilter enabled");
            search_queries.resize(1);
            search_queries[0] = queries;
        }
        else
        {
            log_debug("Applying IBF prefilter");
            std::tie(bin_size, search_queries) = filter_queries(queries, options);
            log_debug("Bin size:", bin_size);
            log_debug("Bucket count:", search_queries.size());
            // auto bucket_sizes = bucket_list | std::views::transform([](auto const & bucket) {
            //     return std::ranges::size(bucket);
            // });
            // auto non_empty_buckets = bucket_sizes | std::views::filter([](size_t const s) { return s > 0;});
            // log_debug("Non-empty bucket count:", std::ranges::distance(non_empty_buckets.begin(), non_empty_buckets.end()));
            // size_t candidates = std::accumulate(bucket_sizes.begin(), bucket_sizes.end(), 0);
            // log_debug("Candidate count:", candidates);
            // log_debug("Bucket sizes:", bucket_sizes);
        }
        end = std::chrono::high_resolution_clock::now();
        log_info("Filter time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        // #TODO: add later
        // std::cout << "bin_size = " << bin_size << "\n";
        // start = std::chrono::high_resolution_clock::now();
        // partitioned_jst_t pjst{std::addressof(jst), bin_size}; // default initialisation
        // std::cout << "bin count = " << pjst.bin_count() << "\n";
        // std::cout << "Run search\n";

        // * filter step with ibf -> {bin_id, {ref_view(query_l)[, ref_view(query_r)], global_query_id}[]}
        // list of {bin_id:queries}
        // partioned_jst[bin_id] -> traverser_model:
            // range_agent{traverser_model, } we can construct this from the model directly.

        // We need to write more information including the reference sequences and the length of the reference sequences.

        // auto print_breakend = [&] (size_t const index) {
        //     auto breakend = *(rcs_store.variants().begin() + index);
        //     log_info("#################");
        //     log_info("breakend index: ", index);
        //     log_info("low breakend: ", libjst::low_breakend(breakend));
        //     log_info("high breakend: ", libjst::high_breakend(breakend));
        //     log_info("#################");
        // };

        // print_breakend(194111);

        start = std::chrono::high_resolution_clock::now();
        // auto query_matches = queries
        //                    | std::views::transform([] (search_query & query) {
        //                         return all_matches{std::move(query)};
        //                    })
        //                    | seqan3::ranges::to<std::vector>();


        using match_positions_t = std::vector<match_position>;
        using bucket_matches_t = std::unordered_map<size_t, match_positions_t>;
        using thread_local_bucket_matches_t = std::vector<bucket_matches_t>;
        thread_local_bucket_matches_t thread_local_matches{};
        thread_local_matches.resize(options.thread_count);

        // now where do we get the chunk size from?
        auto chunked_rcms = rcs_store | libjst::chunk(bin_size);

        size_t bucket_counts{};
        std::ranges::for_each(search_queries, [&] (auto const & bucket) {
            bucket_counts += bucket.size();
        });
        log_info("Total bucket count: ", bucket_counts);

        // parallelising this?
        // go over the bins in dynamic fashion
        #pragma omp parallel for num_threads(options.thread_count) shared(chunked_rcms, thread_local_matches, search_queries, options) schedule(dynamic)
        for (std::ptrdiff_t bin_idx = 0; bin_idx < std::ranges::ssize(chunked_rcms); ++bin_idx)
        { // parallel region
            auto const & bucket_queries = search_queries[bin_idx];
            if (bucket_queries.empty())
                continue;

            bucket_matches_t & local_matches = thread_local_matches[omp_get_thread_num()];
            // Step 1: distribute search:
            log_debug("Local search in bucket: ", bin_idx);
            bucket current_bucket{.base_tree = chunked_rcms[bin_idx],
                                  .needle_list = bucket_queries | std::views::transform([] (search_query const & query) {
                                        return std::views::all(query.value().sequence());
                                   })};
            // std::cout << "std::ranges::size(current_bucket.base_tree.data().source()) = " <<
            //              std::ranges::size(current_bucket.base_tree.data().source()) << "\n";
            // Step 4: apply matching
            log_debug("Initiate searcher");
            bucket_searcher searcher{std::move(current_bucket), options.error_rate};
            // searcher([] () { std::cout << "Found match!\n"; });
            // ATTENTION!!! Not thread safe!
            searcher([&] (std::ptrdiff_t query_idx, match_position position) {
                // log_debug("Record match for query ", query_idx, " at ", position);
                local_matches[bucket_queries[query_idx].key()].push_back(std::move(position));
            });
        }

        end = std::chrono::high_resolution_clock::now();
        log_info("Matching time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        // Step 5: postprocess matches
        start = std::chrono::high_resolution_clock::now();

        // std::vector<search_matches> aligned_matches_list{};
        // aligned_matches_list.reserve(query_matches.size());
        size_t match_count{};
        std::ranges::for_each(thread_local_matches, [&] (bucket_matches_t const & _bucket_matches) {
            for (auto const & match_positions : _bucket_matches) {
                match_count += match_positions.second.size();
            }
        });
            // search_matches aligned_matches{std::move(query_match).query()};
            // match_aligner aligner{rcs_store, aligned_matches.query().value().sequence()};

            // for (match_position pos : query_match.matches()) {
            //     aligned_matches.record_match(aligner(std::move(pos)));
            // }
            // aligned_matches_list.push_back(std::move(aligned_matches));
        std::cout << "match_count: " << match_count << "\n";

        end = std::chrono::high_resolution_clock::now();
        log_info("Aligning time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
        // Filter globally for the best hits.
        // Ater reducing them to the same file sort, uniquify and then erase redundant matches.

        // Step 6: finalise
        start = std::chrono::high_resolution_clock::now();
        // bam_writer writer{rcs_store, options.map_output_file_path};
        // std::ranges::for_each(aligned_matches_list, [&] (search_matches const & matches) {
        //     writer.write_matches(matches);
        // });
        end = std::chrono::high_resolution_clock::now();
        log_info("Writing time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
    }
    catch (std::exception const & ex)
    {
        log_err(ex.what());
        return -1;
    }
    auto global_end = std::chrono::high_resolution_clock::now();
    log_info("Finished mapping [", std::chrono::duration_cast<std::chrono::seconds>(global_end - global_start).count() ,"s]");
    return 0;
}

} // namespace jstmap
