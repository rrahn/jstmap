// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <algorithm>
#include <charconv>
#include <chrono>
#include <ranges>
#include <span>
#include <string>

#include <seqan/vcf_io.h>

#include <seqan3/io/record.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/create/vcf_parser.hpp>
#include <jstmap/create/stripped_vcf_record.hpp>
#include <jstmap/create/serialise_jst.hpp>

namespace jstmap
{
inline reference_t load_base_sequence(std::filesystem::path const & reference_file, std::string_view contig_name)
{
    using namespace std::literals;
    seqan3::sequence_file_input<sequence_input_traits> record_contig_names{reference_file};

    std::string chr_id = /*"chr"s + */std::string{contig_name};
    for (auto && record : record_contig_names) {
        if (record.id().starts_with(chr_id))
            return record.sequence();
    }
    throw std::runtime_error{"Could not find a contig with the name <"s + chr_id + ">!"s};
    return {};
}

inline reference_t reference_for_record(std::filesystem::path const & reference_file_path,
                                        seqan2::VcfFileIn const & vcf_file,
                                        seqan2::VcfRecord const & record) {
    // Load reference sequence!
    return load_base_sequence(reference_file_path,
                              seqan2::toCString(seqan2::contigNames(seqan2::context(vcf_file))[record.rID]));
}

void construct_jst_from_vcf2(std::filesystem::path const & reference_file,
                             std::filesystem::path const & vcf_file_path,
                             std::filesystem::path const & out_file_path)
{
    // Get the application logger.
    auto & log = get_application_logger();

    auto duration = [] (auto const & start) {
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() -
                    start).count();
    };

    // ----------------------------------------------------------------------------
    // Parse the vcf file.
    // ----------------------------------------------------------------------------

    log(verbosity_level::verbose, logging_level::info, "Initialise parsing vcf file ", vcf_file_path);

    // ----------------------------------------------------------------------------
    // Open vcf file handle.
    seqan2::VcfFileIn vcf_file{vcf_file_path.c_str()};
    seqan2::VcfHeader vcf_header{};
    seqan2::readHeader(vcf_header, vcf_file);

    if (seqan2::atEnd(vcf_file))
    {
        log(verbosity_level::standard,
            logging_level::warning,
            "The vcf file ", vcf_file_path, " does not contain any records!");
        return;
    }

    // ----------------------------------------------------------------------------
    // Transform vcf record into intermediate variants.

    auto start = std::chrono::high_resolution_clock::now();

    seqan2::VcfRecord record{};
    seqan2::readRecord(record, vcf_file);
    reference_t reference = reference_for_record(reference_file, vcf_file, record);
    uint32_t haplotype_count = seqan2::length(seqan2::sampleNames(seqan2::context(vcf_file))) * 2;
    std::cout << "haplotype_count: " << haplotype_count << "\n";

    log(verbosity_level::verbose,
        logging_level::info,
        "Loading reference sequence: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Parsing vcf records
    // ----------------------------------------------------------------------------

    start = std::chrono::high_resolution_clock::now();
    rcs_store_t rcs_store{std::move(reference), haplotype_count};
    log_info("Haplotype count: ", haplotype_count);
    size_t record_count{0};
    variant_stat stat{};
    while (!seqan2::atEnd(vcf_file))
    {
        // if (record_count++ % 1000 == 0) {
        //     log_debug("Processing record: ", record_count);
        // }
        stripped_vcf_record tmp_record{vcf_file, rcs_store.variants().coverage_domain(), stat};
        tmp_record.alternatives(rcs_store);
    }

    log_info("Time parsing vcf: ", duration(start), " s");
    log_info("#SNVs: ", stat.snv_count);
    log_info("#InDels: ", stat.indel_count);

    // // ----------------------------------------------------------------------------
    // // Sort the variants.  VCF is sorted by specification?

    // // different model representations -> so the forward_jst is over the reference and the sorted objects model
    // start = std::chrono::high_resolution_clock::now();
    // libjst::variant_store_sorted<variant_store_t> sorted_store{intermediate_variant_store};
    // log(verbosity_level::verbose,
    //     logging_level::info,
    //     "Time sorting jst: ", duration(start), " s");

    // // ----------------------------------------------------------------------------
    // // Resolve ambiguities.
    // start = std::chrono::high_resolution_clock::now();

    // auto first = std::ranges::begin(rcs_store);
    // auto last = std::ranges::end(rcs_store);

    // auto resolve_coverage = [] (auto first, auto last, auto && cmp)
    // {
    //     auto overlap_end = std::ranges::find_if(first, last, std::forward<decltype(cmp)>(cmp));
    //     for (; first != overlap_end; ++first) {
    //         for (auto it = std::ranges::next(first, 1, overlap_end); it != overlap_end; ++it) {
    //                 (libjst::coverage(*it)).and_not(libjst::coverage(*first));
    //         }
    //     }
    //     return first;
    // };

    // auto effective_size = [](auto const &variant) {
    //     return std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
    // };

    // while (first != last) {
    //     // find first adjacent position overlap
    //     first = std::ranges::adjacent_find(first, last, [] (auto const & lhs, auto const & rhs) {
    //         return libjst::position(lhs) == libjst::position(rhs);
    //     });

    //     // Resolve overlaps of insertions at current position overlap
    //     first = resolve_coverage(first, last, [&](auto const & variant) {
    //         return libjst::position(variant) != libjst::position(*first) || effective_size(variant) <= 0;
    //     });

    //     // Resolve overlaps of remaining variants at current position overlap
    //     first = resolve_coverage(first, last, [&](auto const & variant) {
    //         return libjst::position(variant) != libjst::position(*first);
    //     });
    // }

    // log(verbosity_level::verbose,
    //     logging_level::info,
    //     "Time resolving ambiguities in store: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Generate jst model
    // ----------------------------------------------------------------------------

    // start = std::chrono::high_resolution_clock::now();
    // fwd_jst_t fwd_jst{jst};

    // log(verbosity_level::verbose,
    //     logging_level::info,
    //     "Time building forward layer: ", duration(start), " s");

    // ----------------------------------------------------------------------------
    // Serialise rcs_store
    // ----------------------------------------------------------------------------

    start = std::chrono::high_resolution_clock::now();

    serialise(rcs_store, out_file_path);

    log(verbosity_level::verbose,
        logging_level::info,
        "Time serialising jst: ", duration(start), " s");
}

}  // namespace jstmap
