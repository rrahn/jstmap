// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <ranges>
#include <set>
#include <span>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <libjst/journaled_sequence_tree.hpp>

#include "test_utility.hpp" // make_gaped

using namespace std::literals;

using alphabet_t = char;
using shared_event_t = libjst::detail::delta_event_shared<alphabet_t>;
using delta_event_t = typename shared_event_t::delta_event_type;
using substitution_t = typename shared_event_t::substitution_type;
using coverage_t = typename shared_event_t::coverage_type;
using jst_events_t = std::vector<shared_event_t>;

struct traversal_fixture
{
    std::string reference{};
    size_t sequence_count{};
    jst_events_t events{};
    size_t context_size{};
};

struct traversal_test : public ::testing::TestWithParam<traversal_fixture>
{
    using aligned_sequence_t = std::vector<seqan3::gapped<alphabet_t>>;
    using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;
    using context_position_map_t = std::map<std::string_view, std::vector<libjst::context_position>>;

    std::vector<std::string> sequences{}; // The generated sequences from the delta events
    std::vector<alignment_t> alignments{}; // The alignments against the reference sequence generated from the delta events.
    context_position_map_t context_position_map{};

    // Variable to validate a correct traversal.
    int64_t total_context_count{};
    std::vector<libjst::context_position> unknown_locations{};

    void SetUp() override
    {
        generate_alignments();
        generate_context_map();

        // for (std::string seq : sequences)
        //     std::cout << seq << "\n";
    }

    bool all_contexts_enumerated() const
    {
        return total_context_count == 0;
    }

    template <typename position_range_t>
        requires std::same_as<std::ranges::range_value_t<position_range_t>, libjst::context_position>
    bool context_positions_exist(std::string_view context, position_range_t && locations)
    {
        if (std::ranges::empty(locations))
            return true;

        if (auto it = context_position_map.find(context); it != context_position_map.end())
        {
            bool found_all{true};
            for (libjst::context_position const & actual_location : locations)
            {
                size_t erased_elements = std::erase(it->second, actual_location);

                EXPECT_LE(erased_elements, 1u);

                if (erased_elements == 0u)
                {
                    unknown_locations.push_back(actual_location);
                    found_all = false;
                }

                --total_context_count;
            }
            return found_all;
        }
        return  false;
    }

    auto construct_jst() const
    {
        libjst::journaled_sequence_tree<std::string> jst{std::string{GetParam().reference}};

        std::ranges::for_each(alignments, [&] (alignment_t const alignment)
        {
            // seqan3::debug_stream << "Test alignment: " << alignment << "\n";
            jst.add(alignment);
        });

        jst.print_event_queue();
        return jst;
    }

private:

    void generate_alignments()
    {
        // We generate all sequences here from the reference and the events.
        // Then we create a map with position and context.
        sequences.resize(GetParam().sequence_count);
        alignments.resize(GetParam().sequence_count);

        // We construct the sequences from the reference and the given delta map.
        for (unsigned i = 0; i < GetParam().sequence_count; ++i)
        {
            std::string first_sequence = GetParam().reference;
            std::string second_sequence = first_sequence;

            int32_t virtual_offset = 0;
            std::ranges::for_each(GetParam().events, [&] (shared_event_t const & event)
            {
                EXPECT_EQ(event.coverage().size(), GetParam().sequence_count);

                // Continue only if the coverage is true for this sequence.
                if (!event.coverage()[i])
                    return;

                std::visit([&] (auto const & event_kind)
                {
                    size_t const event_position = event.position() + virtual_offset;
                    size_t const insertion_size = event.insertion_size();
                    size_t const deletion_size = event.deletion_size();

                    EXPECT_LE(event_position, first_sequence.size());
                    EXPECT_LE(event_position, second_sequence.size());

                    seqan3::detail::multi_invocable
                    {
                        [&] <typename alphabet_t> (libjst::detail::delta_kind_substitution<alphabet_t> const & s)
                        {
                            // aaaaaaaaa
                            // aaaabbbaa
                            second_sequence.replace(event_position, insertion_size, s.value().data());
                        },
                        [&] <typename alphabet_t> (libjst::detail::delta_kind_insertion<alphabet_t> const & i)
                        {
                            // aaaa--aaaaa
                            // aaaabbaaaaa
                            first_sequence.insert(event_position, insertion_size, '-');
                            second_sequence.insert(second_sequence.begin() + event_position,
                                                   i.value().begin(), i.value().end());
                            virtual_offset += insertion_size;
                        },
                        [&] (libjst::detail::delta_kind_deletion const &)
                        {
                            // aaaaaaaaaaaa
                            // aaaaa----aaa
                            second_sequence.replace(event_position, deletion_size, deletion_size, '-');
                        },
                        [] (...)
                        {
                            FAIL();
                        }
                    } (event_kind);
                }, event.delta_variant());
            });

            // Store the generated alignment.
            alignments[i] = std::pair{libjst::test::make_gapped(first_sequence),
                                      libjst::test::make_gapped(second_sequence)};

            // seqan3::debug_stream << "Alignment " << "i: " << alignments[i] << "\n";

            std::erase(second_sequence, '-');
            sequences[i] = second_sequence;
        }
    }

    void generate_context_map()
    {
        size_t const ctxt_size = GetParam().context_size;

        size_t sequence_index = 0;
        std::ranges::for_each(sequences, [&] (std::string const & sequence)
        {
            std::string_view sv{sequence};
            size_t context_end_position = std::max<int32_t>(sv.size() - ctxt_size + 1, 0);
            assert(context_end_position <= sv.size());
            for (size_t context_position = 0; context_position < context_end_position; ++context_position)
            {
                libjst::context_position context_location{.sequence_id = sequence_index,
                                                          .sequence_position = context_position};
                using value_t = context_position_map_t::value_type;
                value_t insert_value{sv.substr(context_position, ctxt_size), std::vector{context_location}};
                if (auto [it, inserted] = context_position_map.insert(std::move(insert_value)); !inserted)
                    it->second.emplace_back(context_location);

                ++total_context_count;
            }
            ++sequence_index;
        });
    }
};

TEST_P(traversal_test, construct)
{
    auto jst = this->construct_jst();

    EXPECT_EQ(jst.size(), this->sequences.size());

    // TODO: jst-> range of sequence interface.
    // for (size_t i = 0; i < jst.size(); ++i)
    //     EXPECT_RANGE_EQ(jst[i], this->sequences[i]);
}

TEST_P(traversal_test, enumerate_contexts)
{
    auto jst = this->construct_jst();
    auto context_enumerator = jst.context_enumerator_2(GetParam().context_size);
    size_t context_id{};
    for (auto context_it = context_enumerator.begin(); context_it != context_enumerator.end(); ++context_it)
    {
        auto context = *context_it;
        std::string tmp{context.begin(), context.end()};

        std::cout << "context " << context_id++ << ": " << tmp << "\n";
        auto positions = context_it.positions();
        for (auto position : positions)
            std::cout << "\t" << position << "\n";

        EXPECT_TRUE((this->context_positions_exist(tmp, positions))) << "context " << tmp;
    }

    // Verify that all unique contexts have been enumerated and that there is no unknown location.
    EXPECT_TRUE(this->all_contexts_enumerated());

    for (auto && [context, positions] : this->context_position_map)
    {
        if (positions.empty())
            continue;

        std::cout << "Context: " << context;
        for (auto && [id, pos] : positions)
            std::cout << "\t [" << id << ", " << pos << "]";

        std::cout << "\n";
    }

    EXPECT_TRUE(this->unknown_locations.empty());

    for (libjst::context_position & unkown_location : this->unknown_locations)
        std::cout << unkown_location << "\n";
}

// ----------------------------------------------------------------------------
// Test substitutions
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(substitution_1, traversal_test, testing::Values(traversal_fixture
{
    //          0123456
    //               b
    // 0:       aaaa     [0, 0, 0, 0]
    // 1:        aaaa    [1, 1, 1, 1]
    // 2:         aaab   [-, 2, 2, -]
    // 3:          aaba  [-, 3, 3, -]
    // 4:         aaaa   [2, -, -, 2]
    // 5:          aaaa  [3, -, -, 3]
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{5u, substitution_t{"b"s}, coverage_t{0, 1, 1, 0}}
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_event_2, traversal_test, testing::Values(traversal_fixture
{
    //           b
    //          0123456
    // 0        abaa      [0, 0, -, -]
    // 1         baaa     [1, 1, -, -]
    // 2        aaaa      [-, -, 0, 0]
    // 3         aaaa     [-, -, 1, 1]
    // 4          aaaa    [2, 2, 2, 2]
    // 5           aaaa   [3, 3, 3, 3]
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{1u, substitution_t{"b"s}, coverage_t{1, 1, 0, 0}}
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_at_begin, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{0u, substitution_t{"b"s}, coverage_t{1, 1, 0, 0}}
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_at_end, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{6u, substitution_t{"b"s}, coverage_t{1, 0, 0, 1}}
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_at_same_position, traversal_test, testing::Values(traversal_fixture
{
    //seq1      aaabada
    //seq2      aaacaaa
    //seq3      aaabaaa
    //seq4      aaaaaaa
    //             c d

    // 00:      aaab     [0, -, 0, -]
    // 01:       aaba    [1, -, 1, -]
    // 02:        abaa   [2, -, 2, -]
    // 03:         baaa  [3, -, 3, -]
    // 04:      aaac     [-, 0, -, -]
    // 05:       aaca    [-, 1, -, -]
    // 06:        acad   [-, 2, -, -]
    // 07:         cada  [-, 3, -, -]
    // 08:      aaaa     [-, -, -, 0]
    // 09:       aaaa    [-, -, -, 1]
    // 10:        aaad   [-, -, -, 2]
    // 11:         aada  [-, -, -, 3]
    // 12:        aaaa   [-, -, -, -]
    // 13:         aaaa  [-, -, -, -]
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{3u, substitution_t{"b"s}, coverage_t{1, 0, 1, 0}},
        shared_event_t{3u, substitution_t{"c"s}, coverage_t{0, 1, 0, 0}},
        shared_event_t{5u, substitution_t{"d"s}, coverage_t{0, 1, 0, 1}}
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_overlapping, traversal_test, testing::Values(traversal_fixture
{
    //          b c
    //          01234
    //  0:      ba      [ 0, -]
    //  1:      aa      [ -, 0]
    //  2:       ac     [ 1, -]
    //  3:        ca    [ 2, -]
    //  4:       aa     [ -, 1]
    //  5:        aa    [ -, 2]
    //  6:         aa   [ 3, 3]
    .reference{"aaaaa"s},

    .sequence_count{2u},
    .events
    {
        shared_event_t{ 0u, substitution_t{"b"s}, coverage_t{1, 0}},
        shared_event_t{ 2u, substitution_t{"c"s}, coverage_t{1, 0}}
    },
    .context_size{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_overlapping_2, traversal_test, testing::Values(traversal_fixture
{
    //          b  c  d  e  f
    //          0123456789012
    // 00:      baaaa           0: [0, -, -, -, -]
    // 01:      aaaca           0: [-, 0, -, -, -]
    // 02:       aacaa          1: [-, 1, -, -, -]
    // 03:        acaad         2: [-, 2, -, -, -]
    // 04:         caada        3: [-, 3, -, -, -]
    // 05:      aaaaa           0: [-, -, 0, 0, 0]
    // 06:       aaaaa          1: [1, -, 1, 1, 1]
    // 07:        aaaad         2: [-, -, -, -, -]
    // 08:         aaada        3: [-, -, -, -, -]
    // 09:          aadaa       4: [-, 4, -, -, -]
    // 10:           adaaa      5: [-, 5, -, -, -]
    // 11:            daaaa     6: [-, 6, -, -, -]
    // 12:        aaaaa         2: [2, -, 2, 2, 2]
    // 13:         aaaaa        3: [3, -, 3, 3, 3]
    // 14:          aaaaa       4: [4, -, 4, 4, 4]
    // 15:           aaaae      5: [5, -, 5, 5, -]
    // 16:            aaaea     6: [6, -, 6, 6, -]
    // 17:             aaeaa    7: [7, -, 7, 7, -]
    // 18:              aeaaf   8: [-, -, -, 8, -]
    // 19:              aeaaa   8: [8, -, 8, -, -]
    // 20:           aaaaa      5: [-, -, -, -, 5]
    // 21:            aaaaa     6: [-, -, -, -, 6]
    // 22:             aaaaa    7: [-, 7, -, -, 7]
    // 23:              aaaaf   8: [-, 8, -, -, 8]
    // 24:              aaaaa   8: [-, -, -, -, -]
    //          0123456789012
    //                 -----
    //          b  c  d  e  f
    .reference{"aaaaaaaaaaaaa"s},
    .sequence_count{5u},
    .events
    {
        shared_event_t{ 0u, substitution_t{"b"s}, coverage_t{1, 0, 0, 0, 0}},
        shared_event_t{ 3u, substitution_t{"c"s}, coverage_t{0, 1, 0, 0, 0}},
        shared_event_t{ 6u, substitution_t{"d"s}, coverage_t{0, 1, 0, 0, 0}},
        shared_event_t{ 9u, substitution_t{"e"s}, coverage_t{1, 0, 1, 1, 0}},
        shared_event_t{12u, substitution_t{"f"s}, coverage_t{0, 1, 0, 1, 1}}
    },
    .context_size{5u}
}));

INSTANTIATE_TEST_SUITE_P(0_event_and_too_large_context, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events{},
    .context_size{8u}
}));

INSTANTIATE_TEST_SUITE_P(1_substitution_and_too_large_context, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{ 3u, substitution_t{"b"s}, coverage_t{1, 0, 0, 0}},
    },
    .context_size{8u}
}));

INSTANTIATE_TEST_SUITE_P(no_event_and_equal_context_size, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events{},
    .context_size{7u}
}));

INSTANTIATE_TEST_SUITE_P(1_substitution_and_equal_context_size, traversal_test, testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{ 3u, substitution_t{"b"s}, coverage_t{1, 0, 0, 0}},
    },
    .context_size{7u}
}));

INSTANTIATE_TEST_SUITE_P(everything_substituted_and_context_size_4, traversal_test,
testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{1u},
    .events
    {
        shared_event_t{0u, substitution_t{"b"s}, coverage_t{1}},
        shared_event_t{1u, substitution_t{"c"s}, coverage_t{1}},
        shared_event_t{2u, substitution_t{"d"s}, coverage_t{1}},
        shared_event_t{3u, substitution_t{"e"s}, coverage_t{1}},
        shared_event_t{4u, substitution_t{"f"s}, coverage_t{1}},
        shared_event_t{5u, substitution_t{"g"s}, coverage_t{1}},
        shared_event_t{6u, substitution_t{"h"s}, coverage_t{1}},
    },
    .context_size{4u}
}));

INSTANTIATE_TEST_SUITE_P(everything_substituted_and_context_size_1, traversal_test,
testing::Values(traversal_fixture
{
    .reference{"aaaaaaa"s},
    .sequence_count{1u},
    .events
    {
        shared_event_t{0u, substitution_t{"b"s}, coverage_t{1}},
        shared_event_t{1u, substitution_t{"c"s}, coverage_t{1}},
        shared_event_t{2u, substitution_t{"d"s}, coverage_t{1}},
        shared_event_t{3u, substitution_t{"e"s}, coverage_t{1}},
        shared_event_t{4u, substitution_t{"f"s}, coverage_t{1}},
        shared_event_t{5u, substitution_t{"g"s}, coverage_t{1}},
        shared_event_t{6u, substitution_t{"h"s}, coverage_t{1}},
    },
    .context_size{1u}
}));

INSTANTIATE_TEST_SUITE_P(complex_substitutions, traversal_test,
testing::Values(traversal_fixture
{
    .reference{"aaaaaaaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        shared_event_t{0u, substitution_t{"bbbbb"s}, coverage_t{1, 0, 0, 0}},
        shared_event_t{1u, substitution_t{"ccccc"s}, coverage_t{0, 1, 0, 1}},
        shared_event_t{1u, substitution_t{"dd"s}, coverage_t{0, 0, 1, 0}},
        shared_event_t{4u, substitution_t{"cc"s}, coverage_t{0, 0, 1, 0}},
        shared_event_t{6u, substitution_t{"eee"s}, coverage_t{1, 0, 0, 0}},
        shared_event_t{7u, substitution_t{"fff"s}, coverage_t{0, 0, 1, 1}},
        shared_event_t{11u, substitution_t{"g"s}, coverage_t{1, 1, 0, 0}},
    },
    .context_size{1u}
}));
