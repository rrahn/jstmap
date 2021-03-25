// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <iterator>
#include <tuple>

#include <seqan3/range/decorator/gap_decorator.hpp>

#include <jstmap/index/journaled_sequence_tree_builder.hpp>

namespace jstmap
{

template <typename reference_t, typename sequence_t>
auto compress(reference_t const & reference, sequence_t const & sequence)
{
    return std::pair{seqan3::gap_decorator{reference}, seqan3::gap_decorator{sequence}};
}

jst_t build_journaled_sequence_tree(std::vector<raw_sequence_t> && sequences)
{
    assert(!sequences.empty());

    // Move ownership of first sequence to jst.
    jst_t jst{std::move(sequences[0])};

    // Add first sequence as empty alignment against themself.
    jst.add(compress(jst.reference(), jst.reference()));

    // Align remaining sequences against reference sequence and add it to the jst.
    std::for_each(std::next(sequences.begin()), sequences.end(), [&] (auto const & sequence) {
        jst.add(compress(jst.reference(), sequence));
    });

    return jst;
}

}  // namespace jstmap
