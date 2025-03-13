// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libspm/seqan/alphabet.hpp>

#include <libjst/rcms/rcs_store_reversed.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/search/seed_prefix_node_cargo.hpp>

using spm::operator""_dna5;

struct seed_prefix_node_cargo_test : public ::testing::Test
{
    using reverse_rcms_t = libjst::rcs_store_reversed<jstmap::cms_t>;
                               //   * *   *     *
                             // 0 1 2 3 4 5 6 7 8 9
                             // 9 8 7 6 5 4 3 2 1 0
                               //0123456789012345
    jstmap::reference_t source{"AAAACCCCGGGGTTTT"_dna5};
    jstmap::rcs_store_t rcs_store{source, 5};
    reverse_rcms_t reverse_rcs_store{rcs_store.variants()};
    size_t breakends_count = 10;

    virtual void SetUp() override {
        auto dom = rcs_store.variants().coverage_domain();
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{1, 1}, "G"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{3, 1}, "G"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{5, 1}, "T"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{7, 1}, "T"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{9, 1}, "A"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{11, 1}, "A"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{13, 1}, "C"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
        rcs_store.add(jstmap::variant_t{libjst::breakpoint{15, 1}, "C"_dna5, jstmap::coverage_t{{0, 1, 2}, dom}});
    }

    auto make_reverse_tree() const noexcept {
        return reverse_rcs_store | libjst::make_volatile();
    }

    auto make_node(libjst::seek_position position) const {
        auto tree = make_reverse_tree() | libjst::labelled()
                                        | libjst::coloured()
                                        | libjst::merge()
                                        | libjst::seek();
        return tree.seek(position);
    }

    libjst::seek_position to_forward_position(libjst::seek_position reverse_position) const {
        auto reverse_node = make_node(reverse_position);
        auto reverse_tree = make_reverse_tree();
        jstmap::seed_prefix_node_cargo forward_cargo{*reverse_node, reverse_tree};
        return forward_cargo.position();
    }

    template <typename descriptor_t>
    std::string to_string(descriptor_t descriptor) const noexcept {
        std::string tmp{};
        for (bool b : descriptor)
            tmp += (b) ? "1" : "0";

        return tmp;
    }

};

TEST_F(seed_prefix_node_cargo_test, reference_path_at_4) {
    libjst::seek_position reverse_position{};
    reverse_position.reset(4, libjst::breakpoint_end::low);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 4);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            EXPECT_EQ(descriptor, libjst::breakpoint_end::low);
        } else {
            FAIL() << "Expected node from reference path.";
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, reference_path_at_2) {
    libjst::seek_position reverse_position{};
    reverse_position.reset(2, libjst::breakpoint_end::high);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 6);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            EXPECT_EQ(descriptor, libjst::breakpoint_end::high);
        } else {
            FAIL() << "Expected node from reference path.";
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, reference_path_at_0) {
    libjst::seek_position reverse_position{};
    reverse_position.reset(0, libjst::breakpoint_end::low);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 8);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            EXPECT_EQ(descriptor, libjst::breakpoint_end::low);
        } else {
            FAIL() << "Expected node from reference path.";
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, reference_path_at_8) {
    libjst::seek_position reverse_position{};
    reverse_position.reset(8, libjst::breakpoint_end::low);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 0);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            EXPECT_EQ(descriptor, libjst::breakpoint_end::low);
        } else {
            FAIL() << "Expected node from reference path.";
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, alternate_path_at_4) {
    libjst::seek_position reverse_position{};
    reverse_position.initiate_alternate_node(4);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 5);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            FAIL() << "Expected node from alternate path.";
        } else {
            EXPECT_EQ(descriptor.size(), 1);
            EXPECT_EQ(to_string(descriptor), "1");
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, alternate_path_at_4_path_0) {
    libjst::seek_position reverse_position{};
    reverse_position.initiate_alternate_node(4);
    reverse_position.next_alternate_node(false);

    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 5);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            FAIL() << "Expected node from alternate path.";
        } else {
            EXPECT_EQ(descriptor.size(), 1);
            EXPECT_EQ(to_string(descriptor), "1");
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, alternate_path_at_4_path_01) {
    libjst::seek_position reverse_position{};
    reverse_position.initiate_alternate_node(4);
    reverse_position.next_alternate_node(false);
    reverse_position.next_alternate_node(true);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 3);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            FAIL() << "Expected node from alternate path.";
        } else {
            EXPECT_EQ(descriptor.size(), 3);
            EXPECT_EQ(to_string(descriptor), "101");
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, alternate_path_at_4_path_011) {
    libjst::seek_position reverse_position{};
    reverse_position.initiate_alternate_node(4);
    reverse_position.next_alternate_node(false);
    reverse_position.next_alternate_node(true);
    reverse_position.next_alternate_node(true);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 2);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            FAIL() << "Expected node from alternate path.";
        } else {
            EXPECT_EQ(descriptor.size(), 4);
            EXPECT_EQ(to_string(descriptor), "1101");
        }
    });
}

TEST_F(seed_prefix_node_cargo_test, alternate_path_at_4_path_0110) {
    libjst::seek_position reverse_position{};
    reverse_position.initiate_alternate_node(4);
    reverse_position.next_alternate_node(false);
    reverse_position.next_alternate_node(true);
    reverse_position.next_alternate_node(true);
    reverse_position.next_alternate_node(false);
    libjst::seek_position forward_position = to_forward_position(reverse_position);
    EXPECT_EQ(forward_position.get_variant_index(), 2);
    forward_position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
        if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
            FAIL() << "Expected node from alternate path.";
        } else {
            EXPECT_EQ(descriptor.size(), 4);
            EXPECT_EQ(to_string(descriptor), "1101");
        }
    });
}
