// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// #include <seqan3/alphabet/nucleotide/dna4.hpp>
// #include <seqan3/test/expect_range_eq.hpp>
// #include <seqan3/test/performance/sequence_generator.hpp>
// #include <seqan3/core/debug_stream.hpp>

// #include <libcontrib/seqan/alphabet.hpp>

// #include <libjst/container/concept_jst.hpp>
// #include <libjst/container/jst_forward.hpp>
// #include <libjst/container/jst_base.hpp>
// #include <libjst/utility/bit_vector.hpp>
// #include <libjst/variant/variant_snp.hpp>
// #include <libjst/variant/variant_generic.hpp>
// #include <libjst/variant/variant_store_composite.hpp>
// #include <libjst/variant/variant_store_covered.hpp>
// #include <libjst/traversal/lazy_tree.hpp>

// template <typename alphabet_type>
// struct lazy_tree_test : public ::testing::Test
// {
//     using alphabet_t = alphabet_type;
//     using sequence_t = std::vector<alphabet_t>;
//     using snp_variant_t = libjst::snp_variant<alphabet_t>;
//     using generic_variant_t = libjst::generic_variant<alphabet_t>;
//     using coverage_t = libjst::bit_vector<>;

//     using snp_store_t = std::vector<snp_variant_t>;
//     using generic_store_t = std::vector<generic_variant_t>;
//     using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;
//     using covered_store_t = libjst::variant_store_covered<composite_store_t, libjst::bit_vector<>>;

//     using jst_t = libjst::jst_base<sequence_t, covered_store_t>;
//     using fwd_jst_t = libjst::jst_forward<jst_t>;
//     using lazy_tree_t = libjst::lazy_tree<fwd_jst_t>;

//     inline static const std::vector<alphabet_t> base_sequence{seqan3::test::generate_sequence<alphabet_t>(200)};
//     inline static const std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

//     snp_variant_t snp0{4, seqan3::assign_char_to('T', alphabet_t{})};
//     snp_variant_t snp1{44, seqan3::assign_char_to('A', alphabet_t{})};
//     snp_variant_t snp2{112, seqan3::assign_char_to('C', alphabet_t{})};
//     generic_variant_t var0{44, insertion_sequence, 10};
//     generic_variant_t var1{93, insertion_sequence, 0};
//     generic_variant_t var2{154, {}, 1};

//     jst_t jst{this->base_sequence, 4};

//     void SetUp() override
//     {
//         using value_t = std::ranges::range_value_t<covered_store_t>;

//         jst.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}});
//         jst.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}});
//         jst.insert(value_t{this->snp2, coverage_t{1, 0, 0, 1}});
//         jst.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}});
//         jst.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}});
//         jst.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}});
//     }
// };

// using test_types = ::testing::Types<jst::contrib::dna4,
//                                     seqan3::dna4,
//                                     seqan3::dna5
//                                     >;
// TYPED_TEST_SUITE(lazy_tree_test, test_types);

// TYPED_TEST(lazy_tree_test, construction)
// {
//     using fwd_jst_t = typename TestFixture::fwd_jst_t;
//     using lazy_tree_t = typename TestFixture::lazy_tree_t;

//     EXPECT_FALSE(std::is_default_constructible_v<lazy_tree_t>);
//     EXPECT_TRUE(std::is_copy_constructible_v<lazy_tree_t>);
//     EXPECT_TRUE(std::is_nothrow_move_constructible_v<lazy_tree_t>);
//     EXPECT_FALSE(std::is_copy_assignable_v<lazy_tree_t>);
//     EXPECT_FALSE(std::is_move_assignable_v<lazy_tree_t>);
//     EXPECT_TRUE(std::is_destructible_v<lazy_tree_t>);
//     EXPECT_TRUE((std::is_constructible_v<lazy_tree_t, fwd_jst_t const &, size_t>));
// }

// TYPED_TEST(lazy_tree_test, concept)
// {
//     using lazy_tree_t = typename TestFixture::lazy_tree_t;
//     EXPECT_TRUE(std::ranges::input_range<lazy_tree_t>);
// }

// TYPED_TEST(lazy_tree_test, iterate)
// {
//     using fwd_jst_t = typename TestFixture::fwd_jst_t;
//     using lazy_tree_t = typename TestFixture::lazy_tree_t;
//     fwd_jst_t fwd_jst{this->jst};
//     lazy_tree_t tree{fwd_jst, 4};

//     for (auto && node : tree) {
//         seqan3::debug_stream << node.sequence() << "\n";
//     }
// }