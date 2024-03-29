// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <concepts>
#include <algorithm>
#include <stack>
#include <string>



#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/rcms/dna_compressed_multisequence.hpp>
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/rcms/rcs_store_reversed.hpp>

#include "../mock/rcs_store_mock.hpp"

namespace jst::test::left_extended_reverse_tree {

using source_t = std::string;
using variant_t = jst::test::variant<uint32_t, source_t, uint32_t, std::vector<uint32_t>>;

struct fixture {
    source_t source{};
    uint32_t coverage_size{4};
    uint32_t extend_size{};
    std::vector<variant_t> variants{};
    std::vector<source_t> expected_labels{};

    template <typename stream_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, fixture>
    friend stream_t & operator<<(stream_t & stream, this_t &&) {
        stream << "fixture";
        return stream;
    }
};

struct test : public ::testing::TestWithParam<fixture> {
    using coverage_type = libjst::bit_coverage<uint32_t>;
    using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

    using cms_t = libjst::dna_compressed_multisequence<source_t, coverage_type>;
    using cms_value_t = std::ranges::range_value_t<cms_t>;
    using rcs_store_t = libjst::rcs_store<source_t, cms_t>;
    using rcs_store_reverse_t = libjst::rcs_store_reversed<cms_t>;

    rcs_store_t _mock;
    std::unique_ptr<rcs_store_reverse_t> _reversed_mock{};

    void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        coverage_domain_type domain = _mock.variants().coverage_domain();

        std::ranges::for_each(GetParam().variants, [&] (auto var) {
            _mock.add(cms_value_t{libjst::breakpoint{var.position, var.deletion},
                                  var.insertion,
                                  coverage_type{var.coverage, domain}});
        });
        _reversed_mock = std::make_unique<rcs_store_reverse_t>(_mock.variants());
    }

    rcs_store_reverse_t const & get_mock() const noexcept{
        return *_reversed_mock;
    }
};

} // namespace jst::test::left_extended_reverse_tree

using namespace std::literals;

using fixture = jst::test::left_extended_reverse_tree::fixture;
using variant_t = jst::test::left_extended_reverse_tree::variant_t;
struct left_extended_reverse_tree : public jst::test::left_extended_reverse_tree::test
{
    using typename jst::test::left_extended_reverse_tree::test::rcs_store_t;
    using jst::test::left_extended_reverse_tree::test::get_mock;
    using jst::test::left_extended_reverse_tree::test::GetParam;

    auto make_tree() const noexcept {
        auto const & rcs_mock = get_mock();
        return libjst::volatile_tree{rcs_mock} | libjst::labelled() | libjst::left_extend(GetParam().extend_size);
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(left_extended_reverse_tree, root_sink) {
    auto tree = make_tree();

    using tree_t = decltype(tree);

    using node_t = libjst::tree_node_t<tree_t>;
    using cargo_t = libjst::tree_label_t<tree_t>;

    auto to_string = [] (auto seq) -> std::string {
        std::string str;
        for (char c : seq)
            str.push_back(c);
        return str;
    };

    node_t r = libjst::root(tree);

    std::vector<std::string> actual_labels{};
    std::stack<node_t> path{};
    path.push(r);

    std::cout << "Labels: ";
    while (!path.empty()) {
        node_t p = std::move(path.top());
        path.pop();
        cargo_t label = *p;
        std::cout << to_string(label.sequence()) << " " << std::flush;
        actual_labels.push_back(to_string(label.sequence()));

        if (auto c_ref = p.next_ref(); c_ref.has_value()) {
            path.push(std::move(*c_ref));
        }
        if (auto c_alt = p.next_alt(); c_alt.has_value()) {
            path.push(std::move(*c_alt));
        }
    }
    std::cout << "\n";

    std::ptrdiff_t expected_count = std::ranges::ssize(GetParam().expected_labels);
    std::ptrdiff_t actual_count = std::ranges::ssize(actual_labels);
    EXPECT_EQ(expected_count, actual_count);
    for (std::ptrdiff_t i = 0; i < std::min(expected_count, actual_count); ++i)
        EXPECT_EQ(to_string(GetParam().expected_labels[i]), actual_labels[i]) << i;
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------
using namespace std::literals;

INSTANTIATE_TEST_SUITE_P(no_variant, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{},
    .expected_labels{"GGGGAAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv0, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{0}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{"GGGGAAA"s, "AAAC"s, "AAC"s, "AAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv7, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{7}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{""s, "C"s, "CGGGAAAA"s, "GGGGAAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv4, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{4}, .insertion{"C"s}, .deletion{1}, .coverage{0}}
    },
    .expected_labels{"GGG"s, "GGGC"s, "GGCAAAA"s, "GGGGAAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv6, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{4}, .insertion{"C"s}, .deletion{1}, .coverage{0}},
        variant_t{.position{6}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_labels{"G"s, "GT"s, "GTG"s, "GTGC"s, "TGCAAAA"s,
                                                      "GTGGAAAA"s,
                                  "GGG"s, "GGGC"s, "GGCAAAA"s,
                                             "GGGGAAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv5, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{4}, .insertion{"C"s}, .deletion{1}, .coverage{0}},
        variant_t{.position{5}, .insertion{"T"s}, .deletion{1}, .coverage{0, 2}}
    },
    .expected_labels{"GG"s, "GGT"s, "GGT"s, "GGTC"s, "GTCAAAA"s,
                                                        "GGTGAAAA"s,
                                "GGG"s, "GGGC"s, "GGCAAAA"s,
                                             "GGGGAAAA"s}
}));

INSTANTIATE_TEST_SUITE_P(snv4_snv4, left_extended_reverse_tree, testing::Values(fixture{
    .source{"AAAAGGGG"s},
    .extend_size{3},
    .variants{
        variant_t{.position{4}, .insertion{"C"s}, .deletion{1}, .coverage{0}},
        variant_t{.position{4}, .insertion{"T"s}, .deletion{1}, .coverage{1, 2}}
    },
    .expected_labels{"GGG"s, "GGGT"s, "GGTAAAA"s,
                                  "GGG"s, "GGGC"s, "GGCAAAA"s,
                                  "GGGGAAAA"s}
}));
