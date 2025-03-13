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
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include "fixture_base_ibf.hpp"

namespace just::bench
{
    template <typename matcher_t>
    class state_manager {
    private:

        using state_t = spm::matcher_state_t<matcher_t>;
        using state_stack_t = std::stack<state_t>;

        matcher_t _matcher{};
        state_stack_t _states{};

    public:

        constexpr explicit state_manager(matcher_t matcher) noexcept :
            _matcher{std::move(matcher)},
            _states{}
        {}

        constexpr void notify_push() {
            _states.push(_matcher.capture());
        }

        constexpr void notify_pop() {
            assert(!_states.empty());
            _matcher.restore(_states.top());
            _states.pop();
        }
    };

    template <typename capture_t>
    class fixture_resumable_pattern_ibf : public fixture_base_ibf<capture_t>
    {
    private:
        using base_t = fixture_base_ibf<capture_t>;
    public:

        fixture_resumable_pattern_ibf() = default;
        virtual ~fixture_resumable_pattern_ibf() = default;

        template <typename matcher_t>
        void run(::benchmark::State & state, matcher_t matcher)
        {
            auto tree_closure = [] (size_t window_size) {
                return libjst::labelled() | libjst::coloured()
                                          | libjst::trim(window_size - 1)
                                          | libjst::prune()
                                          | libjst::merge();
            };

            auto make_pattern = [&] (auto const &) {  return matcher; };

            state_manager manager{matcher};
            base_t::run(state, make_pattern, tree_closure, [manager] (auto const & tree) mutable {
                libjst::tree_traverser_base path{tree};
                path.subscribe(manager);
                return path;
            });
            this->processed_bytes = base_t::total_bytes(spm::window_size(matcher));
        }
    };
}  // namespace just::bench

