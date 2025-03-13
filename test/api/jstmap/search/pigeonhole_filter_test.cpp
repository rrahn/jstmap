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

struct pigeonhole_filter_test : public ::testing::Test
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
        ++it; // skip first.
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

// TEST_F(pigeonhole_filter_test, source_only) {

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

//     jstmap::bucket test_bucket{.base_tree = flat_tree, .needle_list = reads};
//     jstmap::pigeonhole_filter filter{test_bucket, 0.0};
//     filter([&] (auto &&, auto && finder, auto && needle_position) {
//         auto [index, offset, hit_size] = needle_position;
//         // std::cout << "needle index = " << index << "\n";
//         // std::cout << "needle offset = " << offset << "\n";
//         // std::cout << "needle hit_size = " << hit_size << "\n";
//         // std::cout << "beginPosition(finder) = " << beginPosition(finder) << "\n";
//         // std::cout << "endPosition(finder) = " << endPosition(finder) << "\n";
//         uint32_t seed_size = endPosition(finder) - beginPosition(finder);
//         EXPECT_EQ(seed_size, hit_size);
//         jstmap::reference_t match_sequence{rcs_store.source().begin() + beginPosition(finder),
//                                            rcs_store.source().begin() + endPosition(finder)};
//         // std::cout << "seed_size: " << seed_size << "\n";
//         // std::cout << "match_sequence: " << std::string{match_sequence.begin(), match_sequence.end()} << "\n";
//         auto const & needle = reads[index];
//         jstmap::reference_t seed{std::ranges::next(needle.begin(), offset),
//                                  std::ranges::next(needle.begin(), offset + hit_size)};
//         EXPECT_EQ(match_sequence, seed);
//     });
// }

TEST_F(pigeonhole_filter_test, complete_tree) {
    auto base_tree = rcs_store | libjst::make_volatile();

    // now we need to generate positions and extract reads from this with errors?
    auto [sampled_positions, reads] = generate_reads(base_tree, 100);
    // constexpr std::ptrdiff_t offset{30'000'000};
    // constexpr std::ptrdiff_t read_length{100};

    // std::vector<jstmap::reference_t> reads{};
    // reads.emplace_back(rcs_store.source().begin() + offset,
    //                    rcs_store.source().begin() + offset + read_length);

    jstmap::bucket test_bucket{.base_tree = base_tree, .needle_list = reads};
    jstmap::pigeonhole_filter filter{test_bucket, 0.0};
    filter([&] (auto && cargo, auto && finder, auto && needle_position) {
        auto [index, offset, hit_size] = needle_position;
        // std::cout << "needle index = " << index << "\n";
        // std::cout << "needle offset = " << offset << "\n";
        // std::cout << "needle hit_size = " << hit_size << "\n";
        // std::cout << "beginPosition(finder) = " << beginPosition(finder) << "\n";
        // std::cout << "endPosition(finder) = " << endPosition(finder) << "\n";
        uint32_t seed_size = endPosition(finder) - beginPosition(finder);
        EXPECT_EQ(seed_size, hit_size);
        auto match_start_it = hostIterator(hostIterator(begin(finder))); // iterator to local haystack segment of current node
        std::ptrdiff_t match_start = std::ranges::distance(cargo.path_sequence().begin(), match_start_it);
        jstmap::reference_t match_sequence{cargo.path_sequence().begin() + match_start,
                                           cargo.path_sequence().begin() + match_start + hit_size};
        // std::cout << "seed_size: " << seed_size << "\n";
        // std::cout << "match_sequence: " << std::string{match_sequence.begin(), match_sequence.end()} << "\n";
        auto const & needle = reads[index];
        jstmap::reference_t seed{std::ranges::next(needle.begin(), offset),
                                 std::ranges::next(needle.begin(), offset + hit_size)};
        EXPECT_EQ(match_sequence, seed);
    });
}

// Identität:
// Experte für, Lehrer für
// In Schritte aufteilen
// Format vorgeben
