// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/match_position.hpp>

#include <libspm/matcher/horspool_matcher.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_unsupported.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

TEST(chunked_tree_test, recover_all_labels) {
    using spm::operator""_dna5;
    auto seq = "CACACACTCAGCATCACACAGGTGAACGTGCTGCAGATGCAGGCAGTCTGGCCTCACTGGCTGCCTCCCTCTACCCAGGCTGCCTCCCTGTACCCAGGCT"_dna5;
    spm::horspool_matcher matcher{seq};

    jstmap::rcs_store_t rcsdb = jstmap::load_jst("/Users/rmaerker/Development/jstmap/build/data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst");

    size_t const window_size = spm::window_size(matcher) - 1;
    auto base_tree = libjst::make_volatile(rcsdb);
    auto chunked_base_tree = base_tree | libjst::chunk(1000u, window_size);

    auto tree_adaptor = libjst::labelled()
                      | libjst::coloured()
                      | libjst::trim(window_size)
                      | libjst::prune_unsupported()
                      | libjst::left_extend(window_size)
                      | libjst::merge() // make big nodes
                      | libjst::seek();

    std::vector<jstmap::match_position> labels_base{};
    std::vector<jstmap::match_position> labels_chunked{};

    auto enumerate = [&] (std::vector<jstmap::match_position> & matches, auto && jst) {
        libjst::tree_traverser_base oblivious_path{jst};
        for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
            auto label = *it;
            // seqan3::debug_stream << "Label: " << label.sequence() << "\n";
            matcher(label.sequence(), [&] ([[maybe_unused]] auto && hystk_finder) {
                matches.push_back(
                    jstmap::match_position{.tree_position{label.position()},
                                           .label_offset{std::ranges::ssize(label.sequence()) -
                                                         seqan2::endPosition(hystk_finder)}});
            });
        }
    };

    std::cout << "Enumerating labels from base_tree" << std::flush;
    enumerate(labels_base, base_tree | tree_adaptor);
    std::cout << " -- done\n";

    std::cout << "Enumerating labels from partial tree" << std::flush;
    for (size_t partial_tree_idx = 0; partial_tree_idx < std::ranges::size(chunked_base_tree); ++partial_tree_idx) {
        size_t old_cnt = labels_chunked.size();
        enumerate(labels_chunked, chunked_base_tree[partial_tree_idx] | tree_adaptor);
        if (labels_chunked.size() > old_cnt) {
            std::cout << "Bucket: " << partial_tree_idx << ",";
        }
    }
    std::cout << " -- done\n";

    EXPECT_EQ(labels_base.size(), labels_chunked.size());
    EXPECT_GT(labels_base.size(), 0u);
    for (unsigned idx = 0; idx < labels_base.size(); ++idx) {
        EXPECT_EQ(labels_base[idx], labels_chunked[idx]);
        std::cout << "EXPECTED: " << labels_base[idx] << "\n";
        std::cout << "CURRENT:  " << labels_chunked[idx] << "\n";
    }
}

// Two points unclear!
//  - first, why not expanding correctly into the left bucket
//  - second, why not using full window size after bucket border
//                                                                                                            61000
// Pos:        60906                                                                                         |     61005                       61033                                                    61089
//                                                                                                           |     *                           *                                                        *
// Seek label:                             GTAATTCACACACTCAGCATCACACAGGTGAACGTGCTGCAGATGCAGGCAGTCTGGCCTCACTGG|CTGCCTCCCTCTACCCAGGCTGCCTCCCTGTACCCAGGCTGCCTCCCTGTGCCCAGGCTGCCTCCCTGCGTACAGTCCACCATGCCAGG|XCCCGGAGCA
// Needle:                                       CACACACTCAGCATCACACAGGTGAACGTGCTGCAGATGCAGGCAGTCTGGCCTCACTGG|CTGCCTCCCTCTACCCAGGCTGCCTCCCTGTACCCAGGCT
// Ref label:  GCCATCCCACAGAAGAGAAAACAGACCAGTAATTCACACACTCAGCATCACACAGGTGAACGTGCTGCAGATGCAGGCAGTCTGGCCTCACTGG|CTGCCGCCCTCTACCCAGGCTGCCTCCCTGTACACAGGCTGCCTCCCTGTGCCCAGGCTGCCTCCCTGCGTACAGTCCACCATGCCAGG|GCCCGGAGCATGGTGG

TEST(chunked_tree_test, unwrap_bin_border) {
    using spm::operator""_dna5;
    auto seq = "CACACACTCAGCATCACACAGGTGAACGTGCTGCAGATGCAGGCAGTCTGGCCTCACTGGCTGCCTCCCTCTACCCAGGCTGCCTCCCTGTACCCAGGCT"_dna5;
    spm::horspool_matcher matcher{seq};

    jstmap::rcs_store_t rcsdb = jstmap::load_jst("/Users/rmaerker/Development/jstmap/build/data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.jst");

    // auto chunked_base_tree = libjst::make_volatile(rcsdb) | libjst::chunk(1000u);
    size_t const window_size = spm::window_size(matcher) - 1;
    auto tree_adaptor = libjst::labelled()
                      | libjst::coloured()
                      | libjst::trim(window_size)
                      | libjst::prune_unsupported()
                      | libjst::left_extend(window_size)
                      | libjst::merge() // make big nodes
                      | libjst::seek();

    // either occurrence starts before or after the break position
    // in first case it must be included inside the smaller bucket
    // in second case it must be included inside larger bucket
    auto base_tree = libjst::make_volatile(rcsdb);
    auto search_tree = base_tree | tree_adaptor;
    constexpr size_t variant_index = 772736;

    for (size_t idx = variant_index - 5; idx < variant_index + 4; ++idx) {
        std::cout << "libjst::position(variant["<< idx << "]): " << libjst::position(*std::ranges::next(base_tree.data().variants().begin(), idx)) << "\n";
    }

    auto init_variant = *std::ranges::next(base_tree.data().variants().begin(), variant_index);
    auto snd_variant = *std::ranges::next(base_tree.data().variants().begin(), variant_index + 3);
    auto nxt_variant = *std::ranges::next(base_tree.data().variants().begin(), variant_index + 4);
    std::cout << "libjst::position(init_variant): " << libjst::position(init_variant) << "\n";
    std::cout << "libjst::position(snd_variant): " << libjst::position(snd_variant) << "\n";
    std::cout << "libjst::position(nxt_variant): " << libjst::position(nxt_variant) << "\n";

    libjst::seek_position pos{};
    pos.initiate_alternate_node(variant_index);
    pos.next_alternate_node(false);
    pos.next_alternate_node(false);
    pos.next_alternate_node(true);
    auto target_node = search_tree.seek(pos);
    auto cargo = *target_node;
    auto label = cargo.sequence();
    size_t end_offset = std::ranges::size(label) - 49;
    size_t begin_offset = end_offset - spm::window_size(matcher);
    seqan3::debug_stream << "Seek slice: " << (label | seqan3::views::slice(begin_offset, end_offset)) << "\n";
    seqan3::debug_stream << "Seek label: " << label << "\n";
    auto const & ref = base_tree.data().source();
    size_t ref_begin = libjst::position(init_variant).value() - window_size;
    size_t ref_end = libjst::position(init_variant).value() + 1 + window_size;
    auto ref_slice = ref | seqan3::views::slice(ref_begin, ref_end);
    seqan3::debug_stream << "Ref label:  " << (ref | seqan3::views::slice(ref_begin, ref_end)) << "\n";
    std::cout << "ref_begin: " << ref_begin << "\n";
    std::cout << "ref_end: " << ref_end << "\n";
    std::cout << "std::ranges::size(label): " << std::ranges::size(label) << "\n";
    std::cout << "std::ranges::size(ref_slice): " << std::ranges::size(ref_slice) << "\n";
}

TEST(chunked_tree_test, bin_extension) {
//                                                                                                                   100
// Pos:         0                                                                                                   |     105                         133                                                      189
//                                                                                                                  |     *                           *                                                        *
// Seek label:                                    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA|AAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAA
// Pattern:                                             AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA|AAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAA
// Ref label:   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA|AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    using spm::operator""_dna5;
    auto ref = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAN"_dna5;

    jstmap::rcs_store_t rcsdb{std::move(ref), 4};
    rcsdb.add( 34, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 1, 1, 0});
    rcsdb.add( 55, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 1, 1, 0});
    rcsdb.add( 67, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 1, 1, 0});
    rcsdb.add( 95, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 1, 1, 0});
    rcsdb.add(103, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 1, 1, 0});
    rcsdb.add(105, jstmap::variant_t{'C'_dna5}, jstmap::coverage_t{1, 1, 0, 0});
    rcsdb.add(107, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 0, 1, 1});
    rcsdb.add(129, jstmap::variant_t{'T'_dna5}, jstmap::coverage_t{0, 0, 1, 1});
    rcsdb.add(133, jstmap::variant_t{'G'_dna5}, jstmap::coverage_t{1, 1, 0, 0});
    rcsdb.add(189, jstmap::variant_t{'G'_dna5}, jstmap::coverage_t{1, 1, 0, 0});

    auto needle = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAA"_dna5;
    spm::horspool_matcher matcher{needle};

    auto search = [&] (auto const & jst) {
        std::vector<jstmap::match_position> occurrences{};
        libjst::tree_traverser_base search_path{jst};
        for (auto it = search_path.begin(); it != search_path.end(); ++it) {
            auto label = *it;
            matcher(label.sequence(), [&] ([[maybe_unused]] auto && hystk_finder) {
                occurrences.push_back(
                    jstmap::match_position{.tree_position{label.position()},
                                           .label_offset{std::ranges::ssize(label.sequence()) -
                                                         seqan2::endPosition(hystk_finder)}});
            });
        }
        return occurrences;
    };

    size_t const window_size = spm::window_size(matcher) - 1;
    auto chunked_base_tree = libjst::make_volatile(rcsdb) | libjst::chunk(100u, window_size);
    auto tree_adaptor = libjst::labelled()
                      | libjst::coloured()
                      | libjst::trim(window_size)
                      | libjst::prune_unsupported()
                      | libjst::left_extend(window_size)
                      | libjst::merge() // make big nodes
                      | libjst::seek();

    auto pjst0 = chunked_base_tree[0] | tree_adaptor;
    auto pjst1 = chunked_base_tree[1] | tree_adaptor;

    std::cout << "Search pjst0\n";
    auto occ_jst0 = search(pjst0);
    std::cout << "Search pjst1\n";
    auto occ_jst1 = search(pjst1);

    EXPECT_EQ(std::ranges::size(occ_jst0), 1u);
    EXPECT_EQ(std::ranges::size(occ_jst1), 0u);

    auto verify = [] (auto const & occs) {
        for (auto occ : occs) {
            std::cout << occ << "\n";
        }
    };

    std::cout << "occ_jst0:\n";
    verify(occ_jst0);
    std::cout << "occ_jst1:\n";
    verify(occ_jst1);
}
