// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the node interface to be used with the lazy tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/detail/multi_invocable.hpp>
#include <seqan3/range/views/type_reduce.hpp>

#include <libjst/container/concept_jst.hpp>
#include <libjst/concept.hpp>
#include <libjst/journal.hpp>

namespace libjst
{
    // Clearly needs some refactoring!
    template <journaled_sequence_tree_c jst_t>
    class jst_node
    {
    private:
        enum branch_kind
        {
            base,
            variant
        };

        using store_t = variant_store_t<jst_t const &>;
        using variant_t = std::ranges::range_value_t<store_t>;
        using position_t = variant_position_t<variant_t>;
        using coverage_t = variant_coverage_t<variant_t>;

        using sequence_t = decltype(libjst::base_sequence(std::declval<jst_t const &>()) | seqan3::views::type_reduce);
        using store_view_t = std::views::all_t<store_t const &>;
        using variant_iterator = std::ranges::iterator_t<store_view_t>;

        // must be compatible with insertion type
        // must be validated by concept first.
        using journal_t = journal<position_t, sequence_t>;
        journal_t _journal{}; // the sequence we may need to represent
        coverage_t _coverage{};
        // later we may add the offset? // or via some external interface?

    public:
        jst_node() = default;
        jst_node(jst_node const &) = default;
        jst_node(jst_node &&) = default;
        jst_node & operator=(jst_node const &) = default;
        jst_node & operator=(jst_node &&) = default;
        explicit jst_node(jst_t const &jst, size_t const window_size) :
            _journal{libjst::base_sequence(jst) | seqan3::views::type_reduce}
            // _store{libjst::variant_store(jst) | std::views::all},
            // _next_variant{std::ranges::begin(_store)},
            // _coverage{coverage_t(libjst::size(jst), true)},
            // _window_size{window_size - 1},
            // _base_size{std::ranges::size(libjst::base_sequence(jst))}
        {
            assert(window_size > 0);

            // if (_next_variant != std::ranges::end(_store))
            // {
            //     auto &&next_variant = *_next_variant;
            //     _next = libjst::position(next_variant);
            //     _last = _next + std::ranges::size(libjst::insertion(next_variant)) + _window_size;
            //     //  _base_size + next_variant.insertion_size() - next_variant.deletion_size());
            // }
            // else
            // {
            //     _next = _base_size;
            //     _last = _next;
            // }
        }

        auto sequence() const noexcept
        {
            return _journal.sequence(); // do not care about resumable state inside of here.
            // auto seq = _journal.sequence();
            // // std::cout << "journal sequence: " << std::string{seq.begin(), seq.end()} << "\n";
            // size_t head_position{};
            // auto it = seq.begin();

            // if constexpr (!is_resumable)
            // {
            //     head_position = _first - std::min(_window_size, _first);
            //     assert(seq.size() >= head_position);
            //     assert(_next >= head_position);
            //     it += head_position;
            // }
            // size_t subrange_size = std::min(std::min(_next, _last), seq.size()) - head_position;
            // return std::ranges::subrange{it, it + subrange_size, subrange_size}; // borrowed range
        }

        bool at_end() const noexcept
        {
            return true; //_next >= _last;
        }

        // could be put into its own strategy class.
        constexpr size_t first_position() const noexcept
        {
            return 0;//_first;
        }

        constexpr size_t next_position() const noexcept
        {
            return 0;//_next;
        }

        constexpr size_t last_position() const noexcept
        {
            return 0;//_last;
        }

        template <typename node_t>
        requires std::same_as<std::remove_reference_t<node_t>, jst_node>
        friend auto bifurcate(node_t &&parent) noexcept
            -> std::pair<std::optional<std::remove_reference_t<node_t>>,
                         std::optional<std::remove_reference_t<node_t>>>
        {
            // here we need to split into the different nodes.
            // this is the function that becomes important in our scenario.
            // ------------------
            // create branch node

            std::optional<jst_node> branch_node{std::nullopt};

            // coverage_t child_coverage = parent._coverage & libjst::coverage(*parent._next_variant);
            // if (child_coverage.any())
            // {
            //     jst_node child{};
            //     child._next_variant = parent._next_variant;
            //     child._store = parent._store;
            //     child._coverage = std::move(child_coverage);
            //     child._base_size = parent._base_size;  // size of journal sequence, because we always need to view the entire range.
            //     child._window_size = parent._window_size;
            //     child._first = parent._next;
            //     child._next = parent._last;
            //     child._last = parent._last;
            //     child._journal = parent._journal;
            //     child._kind = branch_kind::variant;
            //     record_sequence_variant(child, *child._next_variant);
            //     // find first branch candidate:
            //     // First find first variant that is not an insertion including the current variant.
            //     child._next_variant =
            //         std::ranges::find_if(++child._next_variant,
            //                              std::ranges::end(child._store),
            //                              [pivot = libjst::position(*parent._next_variant)](auto &&variant)
            //                              {
            //                                  return !libjst::is_insertion(variant) ||
            //                                         libjst::position(variant) != pivot;
            //                              });
            //     // second: if next variant is not already the next valid we need to search it.
            //     auto last_position = [] (auto &&variant)
            //     {
            //         return libjst::position(variant) + libjst::deletion(variant);
            //     };
            //     // Of course! We have to check again, if the variant is at end.
            //     if (child._next_variant != std::ranges::end(child._store) &&
            //         last_position(*parent._next_variant) > libjst::position(*child._next_variant))
            //     {
            //         child._next_variant =
            //             std::ranges::lower_bound(child._next_variant,
            //                                      std::ranges::end(child._store),
            //                                      last_position(*parent._next_variant),
            //                                      std::less<>{},
            //                                      [] (auto &&variant) { return libjst::position(variant); });
            //     }

            //     if (child._next_variant != std::ranges::end(child._store))
            //     {
            //         child._next = parent._next + std::ranges::size(libjst::insertion(*parent._next_variant)) +
            //                       libjst::position(*child._next_variant) - last_position(*parent._next_variant);
            //     }
            //     branch_node = std::move(child);
            // }

            // ------------------
            // create split node

            std::optional<jst_node> split_node{std::nullopt};
            // parent._first = parent._next;
            // auto const &prev_variant = *parent._next_variant;
            // ++parent._next_variant;
            // if (parent._kind == branch_kind::base)
            // { // update base branch node
            //     parent._next = parent._base_size;
            //     parent._last = parent._base_size;
            //     // also wrong! should be the last node set for the base branch.
            //     if (parent._next_variant != std::ranges::end(parent._store))
            //     {
            //         parent._next = libjst::position(*parent._next_variant);
            //         parent._last = parent._next + std::ranges::size(libjst::insertion(*parent._next_variant)) +
            //                        parent._window_size;
            //         //    , parent._base_size +
            //         //                 (*parent._next_variant).event_handle()->insertion_size() -
            //         //                 (*parent._next_variant).event_handle()->deletion_size());
            //     }
            //     split_node = std::move(parent);
            // }
            // else
            // { // update variant branch node
            //     parent._coverage.and_not(libjst::coverage(prev_variant));
            //     if (parent._coverage.any())
            //     { // is there at least one sequence covering the base branch at this site.
            //         if (parent._next_variant != std::ranges::end(parent._store))
            //         { // might end before.
            //             parent._next += libjst::position(*parent._next_variant) - libjst::position(prev_variant);
            //         }
            //         else
            //         {
            //             parent._next = parent._last; // set to last.
            //         }
            //         // assert(prev_variant.position().offset <= next_position);
            //         // parent._next = std::min(parent._next + next_position - prev_variant.position().offset, parent._last);
            //         split_node = std::move(parent);
            //     }
            // }
            return {std::move(branch_node), std::move(split_node)};
        }

    private:
        template <typename sequence_variant_t>
        static void record_sequence_variant(jst_node &node, sequence_variant_t const &variant) noexcept
        {
            if (libjst::is_insertion(variant)) {
                node._journal.record_insertion(node._first, libjst::insertion(variant));
            } else if (libjst::is_deletion(variant)) {
                node._journal.record_deletion(node._first, libjst::deletion(variant));
            } else {
                assert(libjst::is_replacement(variant));
                node._journal.record_substitution(node._first, libjst::insertion(variant));
            }
        }
    };
} // namespace libjst