// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::naive_search.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan3/range/concept.hpp>

#include <libjst/search/state_manager_concept.hpp>
#include <libjst/search/state_manager_single.hpp>

namespace libjst
{

/*!\brief A naive pattern search implementation.
 *
 * \tparam pattern_t The type of the pattern to search; must model std::ranges::view and seqan3::sequence.
 * \tparam state_manager_t The type of the state manager; must model libjst::search_state_manager.
 *
 * \details
 *
 * Implements a naïve search algorithm by scanning the pattern and comparing it characterwise to each position of the
 * text. If a hit was found it calles the given callback with the current iterator/cursor.
 */
template <seqan3::sequence pattern_t,
          search_state_manager state_manager_t = search_state_manager_single<size_t>>
//!\cond
    requires std::ranges::view<pattern_t> && std::ranges::random_access_range<pattern_t>
//!\endcond
class naive_pattern_searcher
{
private:
    pattern_t _pattern{}; //!\brief The pattern to search for.
    state_manager_t _state_manager{}; //!\brief The state manager to use for the search.
public:

    using state_type = size_t; //!< The type of the state.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr naive_pattern_searcher() = default; //!< Default.
    constexpr naive_pattern_searcher(naive_pattern_searcher const &) = default; //!< Default.
    constexpr naive_pattern_searcher(naive_pattern_searcher &&) = default; //!< Default.
    constexpr naive_pattern_searcher & operator=(naive_pattern_searcher const &) = default; //!< Default.
    constexpr naive_pattern_searcher & operator=(naive_pattern_searcher &&) = default; //!< Default.
    ~naive_pattern_searcher() = default; //!< Default.

    /*!\brief Constructs a new naive pattern searcher.
     *
     * \tparam sequence_t The type of the pattern; must model std::ranges::viewable_range and seqan3::sequence.
     * \tparam other_state_manager_t The type of the state manager to use for the search.
     *
     * \param[in] pattern The pattern to search.
     * \param[in] state_manager The state manager to use. Defaults to libjst::search_state_manager_single
     *
     * \details
     *
     * Initialises the pattern by wrapping it in a std::views::all view.
     */
    template <std::ranges::viewable_range sequence_t,
              search_state_manager other_state_manager_t = search_state_manager_single<state_type>>
    //!\cond
        requires (!std::same_as<std::decay_t<sequence_t>, naive_pattern_searcher>) && seqan3::sequence<sequence_t>
    //!\endcond
    naive_pattern_searcher(sequence_t && pattern, other_state_manager_t && state_manager = {}) :
        _pattern{std::views::all(std::forward<sequence_t>(pattern))},
        _state_manager{std::forward<other_state_manager_t>(state_manager)}
    {}
    //!\}

    /*!\brief Invokes the pattern search on the given input range.
     *
     * \tparam haystack_t The type of the haystack to search; must model std::ranges::input_range.
     * \tparam callback_t The type of the callback to invoke on a hit; must model std::invocable with the iterator
     *                    of the haystack.
     *
     * \param[in] haystack The haystack to search the pattern in.
     * \param[in] on_hit The callback to invoke when a hit is found.
     *
     * \details
     *
     * Iterates over the haystack and checks at every position if the pattern can be found withot any errors.
     * If the pattern was found the callback will be invoked with the iterator pointing to the current position.
     */
    template <std::ranges::input_range haystack_t, std::invocable<std::ranges::iterator_t<haystack_t> &> callback_t>
    constexpr void operator()(haystack_t && haystack, callback_t && on_hit) noexcept
    {
        for (auto it = haystack.begin(); it != haystack.end(); ++it)
        {
            std::ranges::range_value_t<haystack_t> value = *it;
            state_type & current_pos = _state_manager.state();
            if (_pattern[current_pos] == value)
            {
                if (++current_pos == std::ranges::size(_pattern))
                {
                    on_hit(it);
                    current_pos = _pattern[0] == value;
                }
            }
            else
            {
                current_pos = _pattern[0] == value;
            }
        }
    }

    //!\brief Returns a reference to the underlying state manager.
    constexpr state_manager_t & state_manager() noexcept
    {
        return _state_manager;
    }
};

/*!\name Type deduction guides
 * \relates libjst::naive_pattern_searcher
 * \{
 */
//!\brief Deduces the pattern type from the constructor argument.
template <std::ranges::viewable_range sequence_t>
naive_pattern_searcher(sequence_t &&) -> naive_pattern_searcher<std::views::all_t<sequence_t>>;

//!\brief Deduces the pattern and the state manager type from the constructor arguments.
template <std::ranges::viewable_range sequence_t, search_state_manager state_manager_t>
naive_pattern_searcher(sequence_t &&, state_manager_t &&)
    -> naive_pattern_searcher<std::views::all_t<sequence_t>, std::decay_t<state_manager_t>>;
//!\}

}  // namespace libjst
