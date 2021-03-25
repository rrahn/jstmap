// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functions to write the results.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

// #include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include <jstmap/search/write_results.hpp>

namespace jstmap
{

// struct sam_record
// {
//     std::string qname{"*"};
//     int16_t flag{};
//     std::string rname{"*"};
//     int32_t pos{0};
//     int32_t mapq{255};
//     std::string cigar{"*"};
//     std::string rnext{"*"};
//     int32_t pnext{0};
//     int64_t tlen{0};
//     std::string seq{"*"};
//     std::string qual{"*"};

//     auto as_record() const
//     {
//         using record_t = seqan3::record<seqan3::type_list<std::string, int32_t>,
//                                            seqan3::fields<field::id,   field::mapq>>;

//         return record_t{qname, mapq};

//         seqan3::detail::for_each([&] (seqan3::field field_id)
//         {
//             switch (field_id)
//             {
//                 case seqan3::field::id: seqan3::get<seqan3::field::id>(tmp) = qname; break;
//                 // case seqan3::field::flag: seqan3::get<seqan3::field::flag>(tmp) = flag; break;
//                 case seqan3::field::ref_id: seqan3::get<seqan3::field::ref_id>(tmp) = rname; break;
//                 case seqan3::field::ref_offset: seqan3::get<seqan3::field::ref_offset>(tmp) = pos; break;
//                 case seqan3::field::mapq: seqan3::get<seqan3::field::mapq>(tmp) = mapq; break;
//                 // case seqan3::field::alignment:
//                 // case seqan3::field::cigar: seqan3::get<seqan3::field::cigar>(tmp) = cigar; break;
//                 // case seqan3::field::mate: seqan3::get<seqan3::field::mate>(tmp) = std::tie(rnext, pnext, tlen); break;
//                 case seqan3::field::seq: seqan3::get<seqan3::field::seq>(tmp) = seq; break;
//                 case seqan3::field::qual: seqan3::get<seqan3::field::qual>(tmp) = qual; break;
//                 // case seqan3::field::tags:
//                 default: throw std::invalid_argument{"The given field is not supported"};
//             }
//         }, field_ids...);

//         return tmp{};
//     }
// };

// Process:
// a) query -> produce results for multiple reference genomes -> dump all data
//
//
void write_results(std::vector<libjst::context_position> && results,
                   std::filesystem::path const & alignment_map_output_path)
{
    // using reference_ids_t = std::vector<std::string>;

    // reference_ids_t reference_ids{};
    // std::vector<size_t> reference_lengths{};

    // We need to write more information including the reference sequences and the length of the reference sequences.
    seqan3::sam_file_output map_output_file{alignment_map_output_path,
                                            //   reference_ids,
                                            //   reference_lengths,
                                            seqan3::fields<seqan3::field::ref_offset>{}};

    std::ranges::for_each(results, [&] (libjst::context_position const & result)
    {
        // map_output_file.push_back(sam_record{.qname = result.query_name(), // only reference?
        //                                      .rname = result.reference_name(), // only reference?
        //                                      .pos = result.reference_position(),  // copy not expensive
        //                                      .mapq = result.score()}}); // mapping quality not expensive

        map_output_file.emplace_back(result.sequence_position);
    });
}

}  // namespace jstmap
