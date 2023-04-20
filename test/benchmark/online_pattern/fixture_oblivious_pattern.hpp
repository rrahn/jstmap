// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include "fixture_base.hpp"

namespace just::bench
{
    template <typename capture_t>
    class fixture_oblivious_pattern : public fixture_base<capture_t>
    {
    private:
        using base_t = fixture_base<capture_t>;
    public:

        fixture_oblivious_pattern() = default;
        virtual ~fixture_oblivious_pattern() = default;

        template <typename matcher_t>
        void run(::benchmark::State & state, matcher_t matcher)
        {
            auto tree = libjst::make_volatile(this->store()) | libjst::labelled()
                                                             | libjst::coloured()
                                                             | libjst::trim(libjst::window_size(matcher) - 1)
                                                             | libjst::prune()
                                                             | libjst::left_extend(libjst::window_size(matcher) - 1)
                                                             | libjst::merge();

            base_t::run(state, matcher, [tree] () { return libjst::tree_traverser_base{tree}; });
            this->processed_bytes = base_t::total_bytes(tree);
        }
    };
}  // namespace just::bench

