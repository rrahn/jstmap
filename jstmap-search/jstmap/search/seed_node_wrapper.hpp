// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seek position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <jstmap/global/match_position.hpp>

namespace jstmap
{
    template <typename extender_t>
    class extension_state_manager {
    private:

        using best_path_match_t = std::pair<match_position, int32_t>;
        using matcher_state_t = spm::matcher_state_t<extender_t>;
        using state_t = std::pair<matcher_state_t, best_path_match_t>;
        using state_stack_t = std::stack<state_t>;


        extender_t & _extender;
        state_stack_t _states{};

    public:

        constexpr explicit extension_state_manager(extender_t & extender) noexcept :
            _extender{extender},
            _states{}
        {
            _states.emplace(matcher_state_t{},
                            best_path_match_t{match_position{}, std::numeric_limits<int32_t>::lowest()});
        }

        constexpr void notify_push() {
            _states.emplace(_extender.capture(), _states.top().second);
        }

        constexpr void notify_pop() {
            assert(!_states.empty());
            _extender.restore(_states.top().first);
            _states.pop();
        }

        constexpr state_t & top() noexcept {
            return _states.top();
        }

        constexpr state_t const & top() const noexcept {
            return _states.top();
        }

        // constexpr void print_state() const noexcept {
        //     seqan3::debug_stream << "[";
        //     std::ranges::for_each(_extender.capture(), [] (uint64_t word) {
        //         seqan3::debug_stream << std::bitset<64>{word} << ", ";
        //     });
        //     seqan3::debug_stream << "]\n";
        // }
    };
}  // namespace jstmap
