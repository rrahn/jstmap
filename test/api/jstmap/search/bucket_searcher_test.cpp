// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <random>
#include <ranges>

#include <libspm/seqan/alphabet.hpp>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/global/match_position.hpp>
#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/bucket.hpp>
#include <jstmap/search/bucket_searcher.hpp>

using spm::operator""_dna5;

struct sample_position {
    libjst::seek_position position;
    std::ptrdiff_t label_offset;
};

struct bucket_searcher_test : public ::testing::Test
{
    jstmap::rcs_store_t rcs_store{};
    size_t max_read_count{100};
    size_t minor_tick_step = max_read_count * 0.01;
    size_t major_tick_step = max_read_count * 0.1;

    virtual void SetUp() override {
        std::filesystem::path file{DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst"};
        rcs_store = jstmap::load_jst(file);
    }

    template <typename sample_tree_t>
    auto sample_positions(sample_tree_t sample_tree, size_t sample_size) {
        std::mt19937 generator{};
        std::uniform_int_distribution<std::ptrdiff_t> generate_step{1, 10000};

        std::ptrdiff_t next_step = generate_step(generator);
        std::vector<sample_position> sampled_positions{};
        libjst::tree_traverser_base traverser{sample_tree};
        auto it = traverser.begin();
        ++it; // skip first
        for (; it != traverser.end(); ++it) {
            auto cargo = *it;
            auto sample_base = cargo.sequence();
            auto sample_end = sample_base.end() - (sample_size - 1);
            // std::cout << "size: " << std::ranges::size(sample_base) << ", ";
            // std::cout << "last: " << std::ranges::distance(sample_base.begin(), sample_end) << "\n";
            for (auto sample_it = sample_base.begin(); sample_it != sample_end; ++sample_it, --next_step) {

                if (next_step < 0) { // now we sample from the current position!
                    if (sampled_positions.size() % minor_tick_step == 0) {
                        if (sampled_positions.size() % major_tick_step == 0) {
                            std::cout << ':' << std::flush;
                        } else {
                            std::cout << '.' << std::flush;
                        }
                    }
                    sampled_positions.emplace_back(cargo.position(), sample_it - sample_base.begin());

                    if (sampled_positions.size() == max_read_count)
                        return sampled_positions;

                    next_step = generate_step(generator);
                }
            }
        }
        std::cout << "\n";
        return sampled_positions;
    }

    template <typename tree_t>
    auto generate_reads(tree_t base_tree, size_t sample_size) {
        auto sample_tree = base_tree | libjst::labelled()
                                     | libjst::coloured()
                                     | libjst::trim(sample_size - 1)
                                     | libjst::prune()
                                     | libjst::left_extend(sample_size - 1)
                                     | libjst::merge()
                                     | libjst::seek();

        auto sampled_positions = sample_positions(sample_tree, sample_size);
        std::vector<jstmap::reference_t> reads{};
        reads.reserve(sampled_positions.size());
        for (auto [position, offset] : sampled_positions) {
            auto node = sample_tree.seek(position);
            auto cargo = *node;
            auto sequence = cargo.sequence();
            jstmap::reference_t sample{sequence.begin() + offset, sequence.begin() + (offset + sample_size)};
            if (std::ranges::any_of(sample, [] (auto symbol) { return static_cast<char>(symbol) == 'N'; }))
                continue;

            reads.push_back(std::move(sample));
        }
        std::cout << "Number of reads " << reads.size() << "\n";
        return std::pair{std::move(sampled_positions), std::move(reads)};
    }

};

// TEST_F(bucket_searcher_test, source_only) {

//     jstmap::reference_t tmp_source{rcs_store.source().begin(), rcs_store.source().end()};
//     jstmap::rcs_store_t rcs_store_empty{std::move(tmp_source), rcs_store.size()};
//     auto flat_tree = rcs_store_empty | libjst::make_volatile();

//     // now we need to generate positions and extract reads from this with errors?
//     auto [sampled_positions, reads] = generate_reads(flat_tree, 100);
//     // constexpr std::ptrdiff_t offset{30'000'000};
//     // constexpr std::ptrdiff_t read_length{100};

//     // std::vector<jstmap::reference_t> reads{};
//     // reads.emplace_back(rcs_store.source().begin() + offset,
//     //                    rcs_store.source().begin() + offset + read_length);

//     auto verify_tree = flat_tree | libjst::labelled()
//                                  | libjst::trim(99u)
//                                  | libjst::merge()
//                                  | libjst::seek();

//     jstmap::bucket test_bucket{.base_tree = flat_tree, .needle_list = reads};
//     jstmap::bucket_searcher searcher{test_bucket, 0.0};
//     searcher([&] (std::ptrdiff_t const needle_index, jstmap::match_position && match) {
//         auto node = verify_tree.seek(match.tree_position);
//         auto cargo = *node;
//         auto root_path_label = cargo.path_sequence();
//         auto match_begin = std::ranges::next(root_path_label.begin(), match.label_offset);
//         auto match_end = std::ranges::next(match_begin, 100);
//         jstmap::reference_t match_segment{match_begin, match_end};
//         EXPECT_EQ(reads[needle_index], match_segment) << "Position: " << match.tree_position
//                                                       << " Offset: " << match.label_offset;
//     });
// }

TEST_F(bucket_searcher_test, complete_tree) {
    auto base_tree = rcs_store | libjst::make_volatile();

    // now we need to generate positions and extract reads from this with errors?
    auto [sampled_positions, tmp_reads] = generate_reads(base_tree, 100);

    std::vector<jstmap::reference_t> reads{tmp_reads};

    auto verify_tree = base_tree | libjst::labelled()
                                 | libjst::merge()
                                 | libjst::seek();

    jstmap::bucket test_bucket{.base_tree = base_tree, .needle_list = reads};
    jstmap::bucket_searcher searcher{test_bucket, 0.0};
    searcher([&] (std::ptrdiff_t const needle_index, jstmap::match_position && match) {
        auto node = verify_tree.seek(match.tree_position);
        auto cargo = *node;
        auto root_path_label = cargo.path_sequence();
        auto match_begin = std::ranges::next(root_path_label.begin(), match.label_offset);
        auto match_end = std::ranges::next(match_begin, 100);
        std::string match_segment{match_begin, match_end};
        std::string needle{reads[needle_index].begin(), reads[needle_index].end()};
        EXPECT_EQ(match_segment, needle) << "Position: " << match.tree_position
                                         << " Offset: " << match.label_offset
                                         << " Needle: " << needle_index;
    });
}
