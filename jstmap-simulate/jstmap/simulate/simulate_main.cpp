// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map simulator.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <algorithm>
#include <chrono>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <libspm/execute/for_each_stream.hpp>
#include <libspm/execute/make_sender.hpp>
#include <libspm/execute/make_stream.hpp>
#include <libspm/execute/run.hpp>
#include <libspm/execute/then.hpp>
#include <libspm/execute/transform_stream.hpp>

#include <jstmap/global/search_matches.hpp>
#include <jstmap/global/all_matches.hpp>
#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/bam_writer.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/search_query.hpp>
#include <jstmap/search/match_aligner.hpp>
#include <jstmap/simulate/options.hpp>
#include <jstmap/simulate/simulate_main.hpp>
#include <jstmap/simulate/read_sampler.hpp>

namespace jstmap
{
template <typename fn_t>
struct debug_fn {

    fn_t _fn;

    debug_fn(fn_t fn) : _fn{fn}
    {
        std::cout << "Initialise debug fn\n";
    }

    debug_fn(debug_fn const & other) : _fn{other._fn}
    {
        std::cout << "Copy debug fn\n";
    }

    debug_fn(debug_fn && other) : _fn{std::move(other._fn)}
    {
        std::cout << "Move debug fn\n";
    }

    debug_fn & operator=(debug_fn const & other)
    {
        std::cout << "Copy assign debug fn\n";
        _fn = other._fn;
        return *this;
    }

    debug_fn & operator=(debug_fn && other) noexcept
    {
        std::cout << "Move assign debug fn\n";
        _fn = std::move(other._fn);
        return *this;
    }

    template <typename ...args_t>
    auto operator()(args_t && ...args) & -> std::invoke_result_t<fn_t &, args_t...>
    {
        return std::invoke(_fn, (args_t&&)args...);
    }

    template <typename ...args_t>
    auto operator()(args_t && ...args) const & -> std::invoke_result_t<fn_t const &, args_t...>
    {
        return std::invoke(_fn, (args_t&&)args...);
    }

    template <typename ...args_t>
    auto operator()(args_t && ...args) && -> std::invoke_result_t<fn_t, args_t...>
    {
        return std::invoke(std::move(_fn), (args_t&&)args...);
    }
};

int simulate_main(seqan3::argument_parser & simulate_parser)
{
    simulate_options options{};

    simulate_parser.add_positional_option(options.input_file,
                                          "The jst to sample reads from.",
                                          seqan3::input_file_validator{{"jst"}});
    simulate_parser.add_positional_option(options.output_file,
                                          "The file containing the sampled reads.",
                                          seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                        {"fa", "fasta"}});

    simulate_parser.add_flag(options.is_quite,
                             'q',
                             "quite",
                             "Disables all logging.",
                             seqan3::option_spec::standard);
    simulate_parser.add_flag(options.is_verbose,
                             'v',
                             "verbose",
                             "Enables expansive debug logging.",
                             seqan3::option_spec::standard);

    simulate_parser.add_option(options.read_size,
                               's',
                               "read-size",
                               "The size of the reads.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{1u, 500u});

    simulate_parser.add_option(options.read_count,
                               'c',
                               "read-count",
                               "The number of reads to sample.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{1u, std::numeric_limits<uint32_t>::max()});

    simulate_parser.add_option(options.error_rate,
                               'e',
                               "error-rate",
                               "The relative error rate.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{0,1});

    try
    {
        simulate_parser.parse();
        if (options.is_quite) {
            get_application_logger().set_verbosity(verbosity_level::quite);
        } else if (options.is_verbose) {
            get_application_logger().set_verbosity(verbosity_level::verbose);
        }

        log_debug("Input file:", options.input_file.string());
        log_debug("Output file:", options.output_file.string());
        log_debug("Read size:", options.read_size);
        log_debug("Read count:", options.read_count);
        log_debug("Error rate:", options.error_rate);

    }
    catch (seqan3::argument_parser_error const & ex)
    {
        log(logging_level::error, ex.what());
        return -1;
    }

    // Load the sequences.
    auto global_start = std::chrono::high_resolution_clock::now();
    try
    {
        log_info("Starting simulation");
        log_debug("Load jst from file", options.input_file);
        auto rcs_store = load_jst(options.input_file);

        log_debug("Initiate simulation");
        // read_sampler sampler{rcs_store};

        auto simulation = execute::make_sender([] (auto const & data) { return read_sampler{data}; }, rcs_store)
                        | execute::then([&] (auto && sampler) {
                            auto start = std::chrono::high_resolution_clock::now();
                            auto sampled_reads = sampler(options.read_count, options.read_size);
                            auto end = std::chrono::high_resolution_clock::now();
                            log_debug("Sampling time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
                            return sampled_reads;
                        })
                        | execute::then([] (auto && sampled_reads) {
                            std::size_t sample_idx{};
                            return execute::make_stream(std::move(sampled_reads))
                                    | execute::transform_stream([sample_idx] (auto && sample_record) mutable -> all_matches {
                                        using namespace std::literals;
                                        auto && [read, match_position] = sample_record;
                                        std::string id_buffer = "jstsim|"s + std::to_string(sample_idx++);
                                        sequence_record_t record{};
                                        record.id() = std::move(id_buffer);
                                        record.sequence() = std::move(read);
                                        all_matches matches{search_query{sample_idx, std::move(record)}};
                                        matches.record_match(std::move(match_position));
                                        return matches;
                                    });
                        })
                        | execute::then([&] (auto && matched_search_query_stream) {
                            return execute::transform_stream(matched_search_query_stream, [&] (auto && sample) {
                                // transfer ownership!
                                search_matches aligned_matches{std::move(sample).query()};
                                match_aligner aligner{rcs_store, aligned_matches.query().value().sequence()};

                                for (match_position pos : sample.matches()) {
                                    aligned_matches.record_match(aligner(std::move(pos)));
                                }
                                return aligned_matches;
                            });
                        })
                        | execute::then([&] (auto && stream) {

                            log_debug("Save sampled reads");

                            auto start = std::chrono::high_resolution_clock::now();
                            seqan3::sequence_file_output sout{options.output_file};
                            bam_writer writer{rcs_store, options.output_file.replace_extension(".sam")};
                            execute::run(execute::for_each_stream((decltype(stream) &&)stream, [&] (search_matches && matches) {
                                sout.push_back(matches.query().value());
                                writer.write_matches(matches);
                            }));

                            auto end = std::chrono::high_resolution_clock::now();
                            log_debug("Saving time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
                        });

        log_debug("Run simulation");
        auto start = std::chrono::high_resolution_clock::now();

        execute::run(simulation);

        auto end = std::chrono::high_resolution_clock::now();
        log_info("Simulation time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
        log_info("Simulation finished successful");
    }
    catch (std::exception const & ex)
    {
        log_err(ex.what());
        return -1;
    }

    return 0;
}

} // namespace jstmap
