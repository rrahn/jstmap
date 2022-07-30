// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "journal_sequence_tree_traversal_test_template.hpp"

#include <libjst/journaled_sequence_tree_forward.hpp>
struct forward_test : public libjst::test::traversal_fixture_base
{};

TEST_P(forward_test, construct)
{
    auto jst = this->construct_jst();

    EXPECT_EQ(jst.size(), this->sequences.size());

    for (size_t i = 0; i < jst.size(); ++i)
        EXPECT_RANGE_EQ(jst.sequence_at(i), this->sequences[i]);
}

struct window_enumerator {

    static constexpr libjst::resume_traversal resume_policy{libjst::resume_traversal::tail_on_breakpoint};
    size_t _window_size;

    template <typename sequence_t, typename callback_t>
    void operator()(sequence_t && sequence, callback_t && callback) {
        size_t sequence_size = std::ranges::size(sequence);
        for (size_t pos = _window_size; pos <= sequence_size; ++pos) {
            callback(sequence | seqan3::views::slice(pos - _window_size, pos));
        }
    }

    size_t window_size() const noexcept {
        return _window_size;
    }
};

struct receiver {

    std::vector<std::string> expected_contexts{};
    size_t count{};

    template <typename sequence_t>
    void set_next(sequence_t && sequence) noexcept {
        // std::cout << std::string{std::ranges::begin(sequence), std::ranges::end(sequence)} << "\n";
        EXPECT_LT(count, expected_contexts.size());
        EXPECT_RANGE_EQ(sequence, expected_contexts[count]);
        if (!std::ranges::equal(sequence, expected_contexts[count])) {
            std::cerr << "count [" << count << "]\n";
        }
        ++count;
    }

    void set_value() noexcept {
        EXPECT_EQ(count, expected_contexts.size());
    }
};

TEST_P(forward_test, enumerate_contexts)
{
    auto jst = this->construct_jst();
    jst.print_event_queue();
    // now we need a suitable test to explore everything.
    libjst::journaled_sequence_tree_forward fwd_jst{std::move(jst)};

    auto sender = fwd_jst.search(window_enumerator{GetParam().context_size});
    auto op = sender.connect(receiver{.expected_contexts = GetParam().expected_contexts});
    op.start();
}

// // ----------------------------------------------------------------------------
// // Test no variants
// // ----------------------------------------------------------------------------

// INSTANTIATE_TEST_SUITE_P(no_variants, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456
//     // 0:       aaaa     [0, 0, 0, 0]
//     // 1:        aaaa    [1, 1, 1, 1]
//     // 2:         aaaa   [2, 2, 2, 2]
//     // 3:          aaaa  [3, 3, 3, 3]
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events{},
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// // ----------------------------------------------------------------------------
// // Test substitutions
// // ----------------------------------------------------------------------------

// INSTANTIATE_TEST_SUITE_P(substitution_1, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456
//     //               b
//     // 0:       aaaa     [0, 0, 0, 0]
//     // 1:        aaaa    [1, 1, 1, 1]
//     // 2:         aaab   [-, 2, 2, -]
//     // 3:          aaba  [-, 3, 3, -]
//     // 4:         aaaa   [2, -, -, 2]
//     // 5:          aaaa  [3, -, -, 3]
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{0, 1, 1, 0}}
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaab", "aaba", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_event_2, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //           b
//     //          0123456
//     // 0        abaa      [0, 0, -, -]
//     // 1         baaa     [1, 1, -, -]
//     // 2        aaaa      [-, -, 0, 0]
//     // 3         aaaa     [-, -, 1, 1]
//     // 4          aaaa    [2, 2, 2, 2]
//     // 5           aaaa   [3, 3, 3, 3]
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 1, 0, 0}}
//     },
//     .context_size{4u},
//     .expected_contexts{"abaa", "baaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 1, 0, 0}}
//     },
//     .context_size{4u},
//     .expected_contexts{"baaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 1}}
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaab", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_at_same_position, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //seq1      aaabada
//     //seq2      aaacaaa
//     //seq3      aaabaaa
//     //seq4      aaaaaaa
//     //             c d

//     // 00:      aaab     [0, -, 0, -]
//     // 01:       aaba    [1, -, 1, -]
//     // 02:        abaa   [2, -, 2, -]
//     // 03:         baaa  [3, -, 3, -]
//     // 04:      aaac     [-, 0, -, -]
//     // 05:       aaca    [-, 1, -, -]
//     // 06:        acad   [-, 2, -, -]
//     // 07:         cada  [-, 3, -, -]
//     // 08:      aaaa     [-, -, -, 0]
//     // 09:       aaaa    [-, -, -, 1]
//     // 10:        aaad   [-, -, -, 2]
//     // 11:         aada  [-, -, -, 3]
//     // 12:        aaaa   [-, -, -, -]
//     // 13:         aaaa  [-, -, -, -]
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"d"s}, libjst::test::coverage_t{0, 1, 0, 1}}
//     },
//     .context_size{4u},
//     .expected_contexts{"aaab", "aaba", "abaa", "baaa", "aaac", "aaca", "acad",
//                        "cada", "aaaa", "aaaa", "aaad", "aada", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_overlapping, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          b c
//     //          01234
//     //  0:      ba      [ 0, -]
//     //  1:      aa      [ -, 0]
//     //  2:       ac     [ 1, -]
//     //  3:        ca    [ 2, -]
//     //  4:       aa     [ -, 1]
//     //  5:        aa    [ -, 2]
//     //  6:         aa   [ 3, 3]
//     .reference{"aaaaa"s},

//     .sequence_count{2u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  0u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  2u}, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{1, 0}}
//     },
//     .context_size{2u},
//     .expected_contexts{"ba", "aa", "ac", "ca", "aa", "aa", "aa"}
// }));

// INSTANTIATE_TEST_SUITE_P(substitution_overlapping_2, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          b  c  d  e  f
//     //          0123456789012
//     // 00:      baaaa           0: [0, -, -, -, -]
//     // 01:      aaaca           0: [-, 0, -, -, -]
//     // 02:       aacaa          1: [-, 1, -, -, -]
//     // 03:        acaad         2: [-, 2, -, -, -]
//     // 04:         caada        3: [-, 3, -, -, -]
//     // 05:      aaaaa           0: [-, -, 0, 0, 0]
//     // 06:       aaaaa          1: [1, -, 1, 1, 1]
//     // 07:        aaaad         2: [-, -, -, -, -]
//     // 08:         aaada        3: [-, -, -, -, -]
//     // 09:          aadaa       4: [-, 4, -, -, -]
//     // 10:           adaaa      5: [-, 5, -, -, -]
//     // 11:            daaaa     6: [-, 6, -, -, -]
//     // 12:        aaaaa         2: [2, -, 2, 2, 2]
//     // 13:         aaaaa        3: [3, -, 3, 3, 3]
//     // 14:          aaaaa       4: [4, -, 4, 4, 4]
//     // 15:           aaaae      5: [5, -, 5, 5, -]
//     // 16:            aaaea     6: [6, -, 6, 6, -]
//     // 17:             aaeaa    7: [7, -, 7, 7, -]
//     // 18:              aeaaf   8: [-, -, -, 8, -]
//     // 19:              aeaaa   8: [8, -, 8, -, -]
//     // 20:           aaaaa      5: [-, -, -, -, 5]
//     // 21:            aaaaa     6: [-, -, -, -, 6]
//     // 22:             aaaaa    7: [-, 7, -, -, 7]
//     // 23:              aaaaf   8: [-, 8, -, -, 8]
//     // 24:              aaaaa   8: [-, -, -, -, -]
//     //          0123456789012
//     //                 -----
//     //          b  c  d  e  f
//     .reference{"aaaaaaaaaaaaa"s},
//     .sequence_count{5u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  0u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  3u}, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  6u}, libjst::test::substitution_t{"d"s}, libjst::test::coverage_t{0, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  9u}, libjst::test::substitution_t{"e"s}, libjst::test::coverage_t{1, 0, 1, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 12u}, libjst::test::substitution_t{"f"s}, libjst::test::coverage_t{0, 1, 0, 1, 1}}
//     },
//     .context_size{5u},
//     .expected_contexts{"baaaa", "aaaca", "aacaa", "acaad", "caada", "aaaaa", "aaaaa", "aaaad", "aaada", "aadaa",
//                        "adaaa", "daaaa", "aaaaa", "aaaaa", "aaaaa", "aaaae", "aaaea", "aaeaa", "aeaaf", "aeaaa",
//                        "aaaaa", "aaaaa", "aaaaa", "aaaaf", "aaaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(0_event_and_too_large_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events{},
//     .context_size{8u},
//     .expected_contexts{}
// }));

// INSTANTIATE_TEST_SUITE_P(1_substitution_and_too_large_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  3u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//     },
//     .context_size{8u},
//     .expected_contexts{}
// }));

// INSTANTIATE_TEST_SUITE_P(no_event_and_equal_context_size, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events{},
//     .context_size{7u},
//     .expected_contexts{"aaaaaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(1_substitution_and_equal_context_size, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset =  3u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//     },
//     .context_size{7u},
//     .expected_contexts{"aaabaaa", "aaaaaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(everything_substituted_and_context_size_4, forward_test,
// testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{1u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::substitution_t{"d"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::substitution_t{"e"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::substitution_t{"f"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"g"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::substitution_t{"h"s}, libjst::test::coverage_t{1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"bcde", "cdef", "defg", "efgh", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(everything_substituted_and_context_size_1, forward_test,
// testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaa"s},
//     .sequence_count{1u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::substitution_t{"d"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::substitution_t{"e"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::substitution_t{"f"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"g"s}, libjst::test::coverage_t{1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::substitution_t{"h"s}, libjst::test::coverage_t{1}},
//     },
//     .context_size{1u},
//     .expected_contexts{"b", "c", "d", "e", "f", "g", "h", "a", "a", "a", "a", "a", "a", "a"}
// }));

// INSTANTIATE_TEST_SUITE_P(complex_substitutions, forward_test,
// testing::Values(
// libjst::test::traversal_fixture
// {
//    //           bbbbb
//    //            ccccc
//    //            dd
//    //               cc
//    //                 eee
//    //                  fff
//    //                      g
//    //           012345678901
//    // 00:       bbbb
//    // 01:        bbbb
//    // 02:         bbba
//    // 03:          bbae
//    // 04:           baee
//    // 05:       accc
//    // 06:        cccc
//    // 07:         cccc
//    // 08:          ccca
//    // 09:           ccaf
//    // 10:            caff
//    // 11:           ccaa
//    // 12:            caaa
//    // 13:       adda
//    // 14:        ddac
//    // 15:         dacc
//    // 16:       aaaa
//    // 17:        aaac
//    // 18:         aacc
//    // 19:          acca
//    // 20:           ccaf
//    // 21:            caff
//    // 22:        aaaa
//    // 23:         aaaa
//    // 24:          aaae
//    // 25:           aaee
//    // 26:            aeee
//    // 27:             eeea
//    // 28:              eeaa
//    // 29:               eaag
//    // 30:          aaaa
//    // 31:           aaaf
//    // 32:            aaff
//    // 33:             afff
//    // 34:              fffa
//    // 35:               ffaa
//    // 36:           aaaa
//    // 37:            aaaa
//    // 38:             aaaa
//    // 39:              aaaa
//    // 40:               aaag
//    // 41:               aaaa
//    //           012345678901
//     .reference{"aaaaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::substitution_t{"bbbbb"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::substitution_t{"ccccc"s}, libjst::test::coverage_t{0, 1, 0, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::substitution_t{"dd"s},    libjst::test::coverage_t{0, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::substitution_t{"cc"s},    libjst::test::coverage_t{0, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::substitution_t{"eee"s},   libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 7u}, libjst::test::substitution_t{"fff"s},   libjst::test::coverage_t{0, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 11u}, libjst::test::substitution_t{"g"s},    libjst::test::coverage_t{1, 1, 0, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"bbbb", "bbbb", "bbba", "bbae", "baee",
//                        "accc", "cccc", "cccc", "ccca", "ccaf", "caff", "ccaa", "caaa",
//                        "adda", "ddac", "dacc",
//                        "aaaa",
//                        "aaac", "aacc", "acca", "ccaf", "caff",
//                        "aaaa", "aaaa",
//                        "aaae", "aaee", "aeee", "eeea", "eeaa", "eaag",
//                        "aaaa",
//                        "aaaf", "aaff", "afff", "fffa", "ffaa",
//                        "aaaa", "aaaa", "aaaa", "aaaa",
//                        "aaag",
//                        "aaaa"}
// }));

// // ----------------------------------------------------------------------------
// // Test insertions
// // ----------------------------------------------------------------------------

// INSTANTIATE_TEST_SUITE_P(single_base_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //
//     //          0123 4567
//     //          aaaa aaaa
//     // 00:      aaaa          [0, 0, 0, 0]
//     // 01:       aaab         [1, 0, 1, 0]
//     // 02:        aaba        [2, 0, 2, 0]
//     // 03:         abaa       [3, 0, 3, 0]
//     // 04:          baaa      [4, 0, 4, 0]
//     // 05:       aaa a        [0, 1, 0, 1]
//     // 06:        aa aa       [0, 2, 0, 2]
//     // 07:         a aaa      [0, 3, 0, 3]
//     // 08:           aaaa     [5, 4, 5, 4]
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaab", "aaba", "abaa", "baaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(single_base_insertion_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //
//     //          01234567
//     //          aaaaaaaa
//     // 00:     baaa         [0, -, -, 0]
//     // 01:      aaaa        [1, 0, 0, 1]
//     // 02:       aaaa       [2, 1, 1, 2]
//     // 03:        aaaa      [3, 2, 2, 3]
//     // 04:         aaaa     [4, 3, 3, 4]
//     // 05:          aaaa    [5, 4, 4, 5]
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"baaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(single_base_insertion_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //
//     //          01234567
//     //          aaaaaaaa
//     // 00:      aaaa          [0, 0, 0, 0]
//     // 01:       aaaa         [1, 1, 1, 1]
//     // 02:        aaaa        [2, 2, 2, 2]
//     // 03:         aaaa       [3, 3, 3, 3]
//     // 04:          aaaa      [4, 4, 4, 4]
//     // 05:           aaab     [5, -, -, 5]
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaab"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_insertions_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          01234567
//     //          aaaaaaaa
//     // 00:      aaaa               [  0,  0,  0,  0]
//     // 01:       aaaa              [  1,  1,  1,  1]
//     // 02:        aaaa             [  2,  2,  2,  2]
//     // 03:         aaaa            [  3,  3,  3,  3]
//     // 04:          aaaa           [  4,  4,  4,  4]
//     // 05:           aaab          [  5,  -,  -,  -]
//     // 06:           aaac          [  -,  5,  -,  -]
//     // 07:            aacc         [  -,  6,  -,  -]
//     // 08:             accc        [  -,  7,  -,  -]
//     // 09:              cccc       [  -,  8,  -,  -]
//     // 10:           aaad          [  -,  -,  5,  -]
//     // 11:            aadd         [  -,  -,  6,  -]
//     // 12:             addd        [  -,  -,  7,  -]
//     // 13:              dddd       [  -,  -,  8,  -]
//     // 14:               dddd      [  -,  -,  9,  -]
//     // 15:                dddd     [  -,  -, 10,  -]
//     // 16:                 dddd    [  -,  -, 11,  -]
//     // 17:                  dddd   [  -,  -, 12,  -]
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"cccc"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"dddddddd"s}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaab", "aaac", "aacc", "accc",
//                        "cccc", "aaad", "aadd", "addd", "dddd", "dddd", "dddd", "dddd", "dddd"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_insertions_overlap, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0   12345678901234567 89
//     //  0:  b___aaddddddddaaaeeea_aagggg
//     //  1:  ccccaaddddddddaaa___a_aa____
//     //  2:  ____aaddddddddaaaeeeafaagggg
//     //  3:  ____aa________aaaeeeafaa____

//     //          01        234   5 67
//     //      ____aa________aaa___a_aa
//     //      b
//     //      cccc
//     //            dddddddd
//     //                       eee
//     //                           f
//     //                              gggg
//     // 00:  b___aadd                      [ 0,  -,  -,  -]
//     // 01:  cccca                         [ -,  0,  -,  -]
//     // 02:   cccaa                        [ -,  1,  -,  -]
//     // 03:    ccaad                       [ -,  2,  -,  -]
//     // 04:     caadd                      [ -,  3,  -,  -]
//     // 05:      aaddd                     [ 1,  4,  0,  -]
//     // 06:       adddd                    [ 2,  5,  1,  -]
//     // 07:        ddddd                   [ 3,  6,  2,  -]
//     // 08:         ddddd                  [ 4,  7,  3,  -]
//     // 09:          ddddd                 [ 5,  8,  4,  -]
//     // 10:           ddddd                [ 6,  9,  5,  -]
//     // 11:            dddda               [ 7, 10,  6,  -]
//     // 12:             dddaa              [ 8, 11,  7,  -]
//     // 13:              ddaaa             [ 9, 12,  8,  -]
//     // 14:               daaae            [10,  -,  9,  -]
//     // 15:               daaa___a         [ -, 13,  -,  -]
//     // 16:      aa________aaa             [ -,  -,  -,  0]
//     // 17:       a________aaae            [ -,  -,  -,  1]
//     // 18:                aaaee           [11,  -, 10,  2]
//     // 19:                 aaeee          [12,  -, 11,  3]
//     // 20:                  aeeea         [13,  -, 12,  4]
//     // 21:                   eeeaf        [ -,  -, 13,  5]
//     // 22:                    eeafa       [ -,  -, 14,  6]
//     // 23:                     eafaa      [ -,  -, 15,  7]
//     // 24:                   eeea_a       [14,  -,  -,  -]
//     // 25:                    eea_aa      [15,  -,  -,  -]
//     // 26:                     ea_aag     [16,  -,  -,  -]
//     // 27:       a________aaa___a         [ -,  -,  -,  -]
//     // 28:                aaa___af        [ -,  -,  -,  -]
//     // 29:                 aa___afa       [ -,  -,  -,  -]
//     // 30:                  a___afaa      [ -,  -,  -,  -]
//     // 31:                      afaag     [ -,  -, 16,  -]
//     // 32:                       faagg    [ -,  -, 17,  -]
//     // 33:                aaa___a_a       [ -, 14,  -,  -]
//     // 34:                 aa___a_aa      [ -, 15,  -,  -]
//     // 35:                  a___a_aag     [ -,  -,  -,  -]
//     // 36:                      a_aagg    [17,  -,  -,  -]
//     // 37:                        aaggg   [18,  -, 18,  -]
//     // 38:                         agggg  [19,  -, 19,  -]

//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"cccc"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::insertion_t{"dddddddd"s}, libjst::test::coverage_t{1, 1, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::insertion_t{"eee"s}, libjst::test::coverage_t{1, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::insertion_t{"f"s}, libjst::test::coverage_t{0, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"gggg"s}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{5u},
//     .expected_contexts{"baadd",
//                        "cccca", "cccaa", "ccaad", "caadd",
//                        "aaddd", "adddd", "ddddd", "ddddd", "ddddd", "ddddd", "dddda", "dddaa", "ddaaa",
//                        "daaae",
//                        "daaaa",
//                        "aaaaa",
//                        "aaaae", "aaaee", "aaeee", "aeeea", "eeeaf", "eeafa", "eafaa",
//                        "eeeaa", "eeaaa", "eaaag",
//                        "aaaaa", "aaaaf", "aaafa", "aafaa", "afaag", "faagg",
//                        "aaaaa", "aaaaa",
//                        "aaaag", "aaagg", "aaggg", "agggg"
//                        }
// }));

// INSTANTIATE_TEST_SUITE_P(insertion_to_get_exactly_one_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //       0 12
//     //      bacaad
//     //  0:  bacaad
//     //  1:  ba_aa_
//     //  2:  _acaa_
//     //  3:  _a_aad
//     //  4:  _a_aa_

//     // 00:  bacaad   [ 0,  -,  -,  -]
//     // 01:  bacaa_   // unsupported
//     // 02:  ba_aad   // unsupported
//     // 03:  ba_aa_   // unsupported
//     // 04:  _acaad   // unsupported
//     // 05:  _acaa_   // unsupported
//     // 06:  _a_aad   // unsupported
//     // 07:  _a_aa_   // unsupported

//     .reference{"aaa"s},
//     .sequence_count{5u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::insertion_t{"c"s}, libjst::test::coverage_t{1, 0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::insertion_t{"d"s}, libjst::test::coverage_t{1, 0, 0, 1, 0}},
//     },
//     .context_size{6u},
//     .expected_contexts{"bacaad"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_insertions_into_empty_reference, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{""s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"b"s}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"cccc"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"dddddddd"s}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"cccc", "dddd", "dddd", "dddd", "dddd", "dddd"}
// }));

// // ----------------------------------------------------------------------------
// // Test deletions
// // ----------------------------------------------------------------------------

// INSTANTIATE_TEST_SUITE_P(single_base_deletion_in_middle, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          aaaaaxaaaa
//     //          aaaaa-aaaa
//     //          0123456789
//     .reference{"aaaaaxbbbb"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaab", "aabb", "abbb",
//                        "aaax", "aaxb", "axbb", "xbbb", "bbbb"}
// }));

// INSTANTIATE_TEST_SUITE_P(single_base_deletion_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"xaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"xaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(single_base_deletion_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaax"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 9u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaax"}
// }));

// INSTANTIATE_TEST_SUITE_P(multi_base_deletion_in_middle, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          aaaaxxxbbb
//     //          aaaa---aaa
//     //          0123456789
//     .reference{"aaaaxxxbbb"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{3}, libjst::test::coverage_t{1, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaab", "aabb", "abbb", "aaax", "aaxx", "axxx", "xxxb", "xxbb", "xbbb"}
// }));

// INSTANTIATE_TEST_SUITE_P(multi_base_deletion_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          ---aaaaaaa
//     //          xxxa
//     //           xxaa
//     //            xaaa
//     //             aaaa
//     //              aaaa
//     //               aaaa
//     //                aaaa
//     //          0123456789
//     .reference{"xxxaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{3}, libjst::test::coverage_t{1, 1, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"xxxa", "xxaa", "xaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(multi_base_deletion_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          aaaaaaa---
//     //          aaaa
//     //           aaaa
//     //            aaaa
//     //             aaaa
//     //              aaax
//     //               aaxx
//     //                axxx
//     //          0123456789
//     .reference{"aaaaaaaxxx"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 7u}, libjst::test::deletion_t{3}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaax", "aaxx", "axxx"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_deletions_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"xyzzaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{0, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"xyzz", "yzza", "zzaa", "zaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_deletions_shortly_after_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     //          a----bbbbb
//     //          ax--zbbbbb
//     //          axy-zbbbbb
//     .reference{"axyyzbbbbb"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{0, 0, 0, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"abbb",
//                        "axzb", "xzbb",
//                        "axyz", "xyzb", "yzbb",
//                        "axyy", "xyyz", "yyzb", "yzbb", "zbbb",
//                        "bbbb", "bbbb"}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_deletions_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"aaaaaaxxyz"s},
//     .sequence_count{6u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{0, 1, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 9u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{0, 0, 0, 1, 0, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaax", "aaxx", "axxy", "xxyz"}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_longer_than_context_in_middle, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"aaaaxxxxbb"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 1}},
//     },
//     .context_size{3u},
//     .expected_contexts{"aaa", "aaa", "aab", "abb",
//                        "aax", "axx", "xxx", "xxx", "xxb", "xbb"}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_longer_than_context_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"xxxxaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 1, 0, 1}},
//     },
//     .context_size{3u},
//     .expected_contexts{"xxx", "xxx", "xxa", "xaa", "aaa", "aaa", "aaa", "aaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_longer_than_context_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"aaaaaaxxxx"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{3u},
//     .expected_contexts{"aaa", "aaa", "aaa", "aaa", "aax", "axx", "xxx", "xxx"}
// }));

// INSTANTIATE_TEST_SUITE_P(one_sequence_deleted, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{10}, libjst::test::coverage_t{1, 0, 0, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(all_sequences_deleted, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{10}, libjst::test::coverage_t{1, 1, 1, 1}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_generating_only_one_context_in_the_middle, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //
//     //      0123456789
//     //      aaaaaaaaaa
//     //  s1: ----aaaa--
//     //  s2: aaaaaaaa--
//     //  s3: ----aaaaaa
//     //  s4: aaaaaaaaaa
//     //
//     // 00:  aaaa          [ -,  0,  -,  0]
//     // 01:   aaaa         [ -,  1,  -,  1]
//     // 02:    aaaa        [ -,  2,  -,  2]
//     // 03:     aaaa       [ -,  3,  -,  3]
//     // 04:      aaaa      [ 0,  4,  0,  4]
//     // 05:       aaaa     [ -,  -,  1,  5]
//     // 06:        aaaa    [ -,  -,  2,  6]
//     .reference{"xxxxaaaazz"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{1, 1, 0, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"xxxx", "xxxa", "xxaa", "xaaa", "aaaa", "aaaz", "aazz"}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_generating_only_one_split_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0123456789
//     //      aabaccaada
//     //  s0: --b-cc--d-
//     //  s1: --b-ccaad-
//     //  s2: --bacc--da
//     //  s3: --baccaada
//     //  s4: aab-cc--da
//     //  s5: aab-ccaad-
//     //  s6: aabacc--d-
//     //  s7: aabaccaada
//     //
//     // 00:  aab-c         [ -, -, -, -, 0, 0, -, -]
//     // 01:   ab-cc        [ -, -, -, -, 1, 1, -, -]
//     // 02:    b-cc--d     [ 0, -, -, -, 2, -, -, -]
//     // 03:    b-cca       [ -, 0, -, -, -, 2, -, -]
//     // 04:  aaba          [ -, -, -, -, -, -, 0, 0]
//     // 05:   abac         [ -, -, -, -, -, -, 1, 1]
//     // 06:    bacc        [ -, -, 0, 0, -, -, 2, 2]
//     // 07:     acc--d     [ -, -, 1, -, -, -, 3, -]
//     // 08:     cc--da     [ -, -, 2, -, 3, -, -, -]
//     // 19:    acca        [ -, -, -, 1, -, -, -, 3]
//     // 10:     ccaa       [ -, 1, -, 2, -, 3, -, 4]
//     // 11:      caad      [ -, 2, -, 3, -, 4, -, 5]
//     // 12:        aada    [ -, -, -, 4, -, -, -, 6]
//     //          0123456789
//     .reference{"aabaccaada"s},
//     .sequence_count{8u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{1, 1, 1, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{1, 0, 1, 0, 1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 9u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0, 0, 1, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aabc", "abcc", "bccd", "bcca",
//                        "aaba", "abac", "bacc",
//                        "accd", "ccda",
//                        "acca", "ccaa", "caad", "aada"},
// }));

// INSTANTIATE_TEST_SUITE_P(larger_deletion_overlaps_smaller_deletions, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0123456789
//     //      aabaccaada
//     //  s0: --b-cc--d-
//     //  s1: --b-ccaad-
//     //  s2: --bacc--da
//     //  s3: --baccaada
//     //  s4: aab-cc--da
//     //  s5: aab-ccaad-
//     //  s6: aa------da
//     //  s7: aa------d-
//     //  s8: aabaccaada

//     // 00:  aa------da    [ -, -, -, -, -, -, 0, -, -]
//     // 01:  aab-c         [ -, -, -, -, 0, 0, -, -, -]
//     // 02:   ab-cc        [ -, -, -, -, 1, 1, -, -, -]
//     // 03:    b-cc--d     [ 0, -, -, -, 2, -, -, -, -]
//     // 04:    b-cca       [ -, 0, -, -, -, 2, -, -, -]
//     // 05:  aaba          [ -, -, -, -, -, -, -, -, 0]
//     // 06:   abac         [ -, -, -, -, -, -, -, -, 1]
//     // 07:    bacc        [ -, -, 0, 0, -, -, -, -, 2]
//     // 08:     acc--d     [ -, -, 1, -, -, -, -, -, -]
//     // 09:     cc--da     [ -, -, 2, -, 3, -, -, -, -]
//     // 10:    acca        [ -, -, -, 1, -, -, -, -, 3]
//     // 11:     ccaa       [ -, 1, -, 2, -, 3, -, -, 4]
//     // 12:      caad      [ -, 2, -, 3, -, 4, -, -, 5]
//     // 13:        aada    [ -, -, -, 4, -, -, -, -, 6]
//     .reference{"aabaccaada"s},
//     .sequence_count{9u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u},
//                                      libjst::test::deletion_t{2},
//                                      libjst::test::coverage_t{1, 1, 1, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u},
//                                      libjst::test::deletion_t{6},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 1, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 3u},
//                                      libjst::test::deletion_t{1},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u},
//                                      libjst::test::deletion_t{2},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 9u},
//                                      libjst::test::deletion_t{1},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"aada", "aabc", "abcc", "bccd", "bcca",
//                        "aaba", "abac", "bacc", "accd", "ccda",
//                        "acca", "ccaa", "caad", "aada"}
// }));

// INSTANTIATE_TEST_SUITE_P(small_deletions_behind_each_other, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0123456789
//     //      baccaaaaaa
//     //  s0: -a--aaaaaa
//     //  s1: -accaaaaaa
//     //  s2: ba--aaaaaa
//     //  s3: baccaaaaaa

//     // 00:  ba--aa       [ -, -, 0, -]
//     // 01:   a--aaa      [ 0, -, 1, -]
//     // 02:  bacc         [ -, -, -, 0]
//     // 02:   acca        [ -, 0, -, 1]
//     // 03:    ccaa       [ -, 1, -, 2]
//     // 04:     caaa      [ -, 2, -, 3]
//     // 05:      aaaa     [ 1, 3, 2, 4]
//     // 06:       aaaa    [ 2, 4, 3, 5]
//     // 06:        aaaa   [ 3, 5, 4, 6]

//     .reference{"baccaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::deletion_t{2}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"baaa", "aaaa", "bacc", "acca", "ccaa", "caaa", "aaaa", "aaaa", "aaaa"}
// }));

// // ----------------------------------------------------------------------------
// // Test mixed variants
// // ----------------------------------------------------------------------------

// INSTANTIATE_TEST_SUITE_P(insertion_at_begin_followed_by_deletion_of_entire_reference, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //           0123456789
//     //      bbbbbaaaaaaaaaa
//     //  s0: bbbbb----------
//     //  s1: bbbbbaaaaaaaaaa
//     //  s2: _____----------
//     //  s3: _____aaaaaaaaaa

//     // 00:  bbbb            [ 0, 0, -, -]
//     // 01:   bbbb           [ 1, 1, -, -]
//     // 02:    bbba          [ -, 2, -, -]
//     // 03:     bbaa         [ -, 3, -, -]
//     // 04:      baaa        [ -, 4, -, -]
//     // 05:       aaaa       [ -, 5, -, 0]
//     // 06:        aaaa      [ -, 6, -, 1]
//     // 07:         aaaa     [ -, 7, -, 2]
//     // 08:          aaaa    [ -, 8, -, 3]
//     // 09:           aaaa   [ -, 9, -, 4]
//     // 10:            aaaa  [ -,10, -, 5]
//     // 11:             aaaa [ -,11, -, 6]

//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"bbbbb"s}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{10}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"bbbb", "bbbb", "bbba", "bbaa", "baaa",
//                        "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(insertion_at_begin_followed_by_deletion_without_valid_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //          0123456789
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"bbb"s}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{10}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"bbba", "bbaa", "baaa",
//                        "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa"}
// }));

// INSTANTIATE_TEST_SUITE_P(insertion_at_begin_followed_by_deletion_with_one_valid_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaax"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"bbb"s}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{9}, libjst::test::coverage_t{1, 0, 1, 0}},
//     },
//     .context_size{4u},
//     .expected_contexts{"bbbx", "bbba", "bbaa", "baaa",
//                        "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaaa", "aaax"}
// }));

INSTANTIATE_TEST_SUITE_P(two_insertions_with_preceding_and_trailing_deletion, forward_test, testing::Values(
libjst::test::traversal_fixture
{
    //          01   234   56
    // s0:      aa---iii---cc
    // s1:      aa---iiibbbcc
    // s2:      aaxxxjjj---cc
    // s3:      aaxxxjjjbbbcc
    // s4:      aa---___---cc
    // s5:      aa---___bbbcc
    // s6:      aaxxx___---cc
    // s7:      aaxxx___bbbcc
    //          0123456789
    .reference{"aaxxxbbbcc"s},
    .sequence_count{8u},
    .events
    {
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::deletion_t{3},       libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::insertion_t{"iii"s}, libjst::test::coverage_t{1, 1, 0, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::insertion_t{"jjj"s}, libjst::test::coverage_t{0, 0, 1, 1, 0, 0, 0, 0}},
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::deletion_t{3},       libjst::test::coverage_t{1, 0, 1, 0, 1, 0, 1, 0}},
    },
    .context_size{4u},

// 0 [b: ([idx: 0, pos: 2], del: 3) ~   <11000100>]
// 1 [b: ([idx: 0, pos: 2], del: 6) ~   <00001000>]
// 2 [b: ([idx: 0, pos: 5], ins: iii) ~ <11000000>]
// 3 [b: ([idx: 0, pos: 5], ins: jjj) ~ <00110000>]
// 4 [b: ([idx: 0, pos: 5], del: 3) ~   <10100010>]
    .expected_contexts{"aaii", "aiii", "aacc", "aabb", "abbb", // first deletion
                       "aaxx", "axxx", "xxxj", "xxjj", "xjjj", "jjjc", "jjcc", "jjjb", "jjbb", "jbbb", // second insertion (first one is not eligible)
                       "xxxc", "xxcc", // second deletion
                       "xxxb", "xxbb", "xbbb", "bbbc", "bbcc"} // remaining base branch
}));

// INSTANTIATE_TEST_SUITE_P(overlapping_insertion_deletion_substitution_at_begin, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{5u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{1, 1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 0, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::substitution_t{"q"s}, libjst::test::coverage_t{0, 1, 1, 0, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(overlapping_insertion_deletion_substitution_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      01234
//     //      aaaaa
//     //  s0: aaaa-i
//     //  s1: aaaaqi
//     //  s2: aaaaq
//     //  s3: aaaa-
//     //  s4: aaaaa

//     // 00:  aaaa       [ 0, 0, 0, 0, 0]
//     // 01:   aaaq      [ -, 1, 1, -, -]
//     // 02:    aaqi     [ -, 2, -, -, -]
//     // 03:   aaa-i     [ 1, -, -, -, -]
//     // 04:   aaaa      [ -, -, -, -, 1]

//     .reference{"aaaaa"s},
//     .sequence_count{5u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 0, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::substitution_t{"q"s}, libjst::test::coverage_t{0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{1, 1, 0, 0, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(deletion_at_end_without_subsequent_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(longer_deletion_at_end_without_subsequent_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(longer_split_deletion_at_end_with_subsequent_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"ii"s}, libjst::test::coverage_t{1, 1, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(longer_split_deletion_at_end_without_subsequent_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u}, libjst::test::deletion_t{1}, libjst::test::coverage_t{1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"ii"s}, libjst::test::coverage_t{0, 0, 0, 1}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(longer_deletion_and_substitution_with_insertion_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"qqq"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{1, 1, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(longer_deletion_and_substitution_without_insertion_at_end, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaa"s},
//     .sequence_count{4u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u}, libjst::test::deletion_t{4}, libjst::test::coverage_t{1, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"qqq"s}, libjst::test::coverage_t{0, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u}, libjst::test::insertion_t{"i"s}, libjst::test::coverage_t{0, 0, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(three_insertions_with_multiple_preceding_and_trailing_events, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{8u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u},
//                                      libjst::test::substitution_t{"pppp"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"ii"s},
//                                      libjst::test::coverage_t{1, 0, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"jjj"s},
//                                      libjst::test::coverage_t{0, 1, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"k"s},
//                                      libjst::test::coverage_t{0, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::substitution_t{"qq"s},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(three_insertions_with_multiple_preceding_and_trailing_events_and_final_insertion, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"aaaaaaaaaa"s},
//     .sequence_count{16u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u},
//                                      libjst::test::substitution_t{"pppp"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 2u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"ii"s},
//                                      libjst::test::coverage_t{1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"jjj"s},
//                                      libjst::test::coverage_t{0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::insertion_t{"k"s},
//                                      libjst::test::coverage_t{0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::substitution_t{"qq"s},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 9u},
//                                      libjst::test::insertion_t{"llll"},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(insertion_in_middle_surrounded_by_deletion_with_one_valid_context, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0123  456789
//     //      xaaabbaaaaay
//     //  s0: x---bb-----y
//     //  s1: x---bbaaaaay
//     //  s2: x---__-----y
//     //  s3: x---__aaaaay
//     //  s4: xaaabb-----y
//     //  s5: xaaabbaaaaay
//     //  s6: xaaa__-----y
//     //  s7: xaaa__aaaaay
//     //                      0  1  2  3  4  5  6  7
//     // 00:  x---bb-----y  [ 0, -, -, -, -, -, -, -]
//     // 01:  x---bba       [ -, 0, -, -, -, -, -, -]
//     // 02:  x---__aaa     [ -, -, -, 0, -, -, -, -]
//     // 03:  xaaa          [ -, -, -, -, 0, 0, 0, 0]
//     // 04:   aaab         [ -, -, -, -, 1, 1, -, -]
//     // 05:    aabb        [ -, -, -, -, 2, 2, -, -]
//     // 06:     abb-----y  [ -, -, -, -, 3, -, -, -]
//     // 07:     abba       [ -, -, -, -, -, 3, -, -]
//     // 08:      bbaa      [ -, 1, -, -, -, 4, -, -]
//     // 09:       baaa     [ -, 2, -, -, -, 5, -, -]
//     // 10:   aaa__-----y  [ -, -, -, -, -, -, 1, -]
//     // 11:   aaa__a       [ -, -, -, -, -, -, -, 1]
//     // 12:    aa__aa      [ -, -, -, -, -, -, -, 2]
//     // 13:     a__aaa     [ -, -, -, -, -, -, -, 3]
//     // 14:        aaaa    [ -, 3, -, 1, -, 6, -, 4]
//     // 15:         aaaa   [ -, 4, -, 2, -, 7, -, 5]
//     // 16:          aaay  [ -, 5, -, 3, -, 8, -, 6]
//     .reference{"xaaaaaaaay"s},
//     .sequence_count{8u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 1, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::insertion_t{"bb"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::deletion_t{5},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 1, 0, 1, 0}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(insertion_at_end_and_begin_of_substitutions_and_deletions, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     //      0123    45    6789
//     //      xaaa____bb____cccy
//     //  s0: x---ii__qq____qqqy
//     //  s1: x---ii__bbkkkkcccy
//     //  s2: x---jjjjqq____qqqy
//     //  s3: x---jjjjbbkkkkcccy
//     //  s4: xaaaii__--____---y
//     //  s5: xaaaii__bb____cccy
//     //  s6: xaaajjjj--____---y
//     //  s7: xaaajjjjbb____ccrr
//     //                             0   1   2   3   4   5   6   7
//     // 00:  x---ii__q           [ 00,  -,  -,  -,  -,  -,  -,  -]
//     // 01:  x---ii__b           [  -, 00,  -,  -,  -,  -,  -,  -]
//     // 02:  x---jjj             [  -,  -, 00, 00,  -,  -,  -,  -]
//     // 03:  xaaa                [  -,  -,  -,  -, 00, 00, 00, 00]
//     // 04:   aaai               [  -,  -,  -,  -, 01, 01,  -,  -]
//     // 05:    aaii              [  -,  -,  -,  -, 02, 02,  -,  -]
//     // 06:     aii__q           [  -,  -,  -,  -,  -,  -,  -,  -] // unsupported branch
//     // 07:      ii__qq          [ 01,  -,  -,  -,  -,  -,  -,  -]
//     // 08:       i__qq____q     [ 02,  -,  -,  -,  -,  -,  -,  -]
//     // 09:     aii__--____---y  [  -,  -,  -,  -, 03,  -,  -,  -]
//     // 10:     aii__b           [  -,  -,  -,  -,  -, 03,  -,  -]
//     // 11:      ii__bb          [  -, 01,  -,  -,  -, 04,  -,  -]
//     // 12:       i__bbk         [  -, 02,  -,  -,  -,  -,  -,  -]
//     // 13:       i__bb____c     [  -,  -,  -,  -,  -, 05,  -,  -]
//     // 14:   aaaj               [  -,  -,  -,  -,  -,  -, 01, 01]
//     // 15:    aajj              [  -,  -,  -,  -,  -,  -, 02, 02]
//     // 16:     ajjj             [  -,  -,  -,  -,  -,  -, 03, 03]
//     // 17:      jjjj            [  -,  -, 01, 01,  -,  -, 04, 04]
//     // 18:       jjjq           [  -,  -, 02,  -,  -,  -,  -,  -]
//     // 19:        jjqq          [  -,  -, 03,  -,  -,  -,  -,  -]
//     // 20:         jqq____q     [  -,  -, 04,  -,  -,  -,  -,  -]
//     // 21:       jjj--____---y  [  -,  -,  -,  -,  -,  -, 05,  -]
//     // 22:       jjjb           [  -,  -,  -, 02,  -,  -,  -, 05]
//     // 23:        jjbb          [  -,  -,  -, 03,  -,  -,  -, 06]
//     // 24:         jbbk         [  -,  -,  -, 04,  -,  -,  -,  -]
//     // 25:         jbb____c     [  -,  -,  -,  -,  -,  -,  -, 07]
//     // 26:   aaa____q           [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported branch
//     // 27:    aa____qq          [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported branch
//     // 28:     a____qq____q     [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported branch
//     // 29:          qq____qq    [ 03,  -, 05,  -,  -,  -,  -,  -]
//     // 30:           q____qqq   [ 04,  -, 06,  -,  -,  -,  -,  -]
//     // 31:                qqqy  [ 05,  -, 07,  -,  -,  -,  -,  -]
//     // 32:   aaa____--____---y  [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported branch
//     // 33:   aaa____b           [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported base
//     // 34:    aa____bb          [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported base
//     // 35:     a____bbk         [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported branch
//     // 36:          bbkk        [  -, 03,  -, 05,  -,  -,  -,  -]
//     // 37:           bkkk       [  -, 04,  -, 06,  -,  -,  -,  -]
//     // 38:            kkkk      [  -, 05,  -, 07,  -,  -,  -,  -]
//     // 39:             kkkc     [  -, 06,  -, 08,  -,  -,  -,  -]
//     // 40:              kkcc    [  -, 07,  -, 09,  -,  -,  -,  -]
//     // 41:               kccc   [  -, 08,  -, 10,  -,  -,  -,  -]
//     // 42:     a____bb____c     [  -,  -,  -,  -,  -,  -,  -,  -]  // unsupported base
//     // 43:          bb____cc    [  -,  -,  -,  -,  -, 06,  -, 08]
//     // 44:           b____ccr   [  -,  -,  -,  -,  -,  -,  -, 09]
//     // 45:                ccrr  [  -,  -,  -,  -,  -,  -,  -, 10]
//     // 46:           b____ccc   [  -,  -,  -,  -,  -, 07,  -,  -]
//     // 47:                cccy  [  -, 09,  -, 11,  -, 08,  -,  -]

//     .reference{"xaaabbcccy"s},
//     .sequence_count{8u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 1, 1, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::insertion_t{"ii"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::insertion_t{"jjjj"s},
//                                      libjst::test::coverage_t{0, 0, 1, 1, 0, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::substitution_t{"qqqqq"s},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::deletion_t{5},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u},
//                                      libjst::test::insertion_t{"kkkk"s},
//                                      libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u},
//                                      libjst::test::substitution_t{"rr"s},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 1}},
//     },
//     .context_size{4u}
// }));

// INSTANTIATE_TEST_SUITE_P(multiple_overlapping_and_nested_variants, forward_test, testing::Values(
// libjst::test::traversal_fixture
// {
//     .reference{"xaaabbcccy"s},
//     .sequence_count{8u},
//     .events
//     {
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u},
//                                      libjst::test::insertion_t{"f"s},
//                                      libjst::test::coverage_t{1, 0, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u},
//                                      libjst::test::insertion_t{"gg"s},
//                                      libjst::test::coverage_t{0, 1, 0, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u},
//                                      libjst::test::insertion_t{"hhh"s},
//                                      libjst::test::coverage_t{0, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 0u},
//                                      libjst::test::substitution_t{"pppp"s},
//                                      libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 1u},
//                                      libjst::test::deletion_t{3},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::insertion_t{"ii"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::insertion_t{"jjjj"s},
//                                      libjst::test::coverage_t{0, 0, 1, 1, 0, 0, 1, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::substitution_t{"qqqqq"s},
//                                      libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 4u},
//                                      libjst::test::deletion_t{5},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 1, 0, 1, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 6u},
//                                      libjst::test::insertion_t{"kkkk"s},
//                                      libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 8u},
//                                      libjst::test::substitution_t{"rr"s},
//                                      libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 1}},
//         libjst::test::shared_event_t{libjst::test::position_t{.offset = 10u},
//                                      libjst::test::insertion_t{"lll"s},
//                                      libjst::test::coverage_t{1, 1, 0, 0, 0, 1, 0, 1}},
//     },
//     .context_size{4u}
// }));