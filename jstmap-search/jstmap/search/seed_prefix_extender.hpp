// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides right extender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <stack>

#include <libspm/matcher/myers_prefix_matcher_restorable.hpp>
#include <libjst/rcms/rcs_store_reversed.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/global/match_position.hpp>
#include <jstmap/search/seed_node_wrapper.hpp>
#include <jstmap/search/seed_extension_tree.hpp>
#include <jstmap/search/seed_prefix_node_cargo.hpp>
#include <jstmap/search/seed_prefix_seek_position.hpp>

namespace jstmap
{
    template <typename base_tree_t, typename needle_t>
    class seed_prefix_extender
    {
        using reverse_rcs_t = decltype(libjst::rcs_store_reversed{std::declval<base_tree_t &&>().data().variants()});
        using tree_t = decltype(std::declval<reverse_rcs_t>() | libjst::make_volatile());
        using reverse_needle_t = decltype(std::declval<needle_t&&>() | std::views::reverse);

        base_tree_t const & _base_tree;
        reverse_rcs_t _reverse_rcms;
        tree_t _reverse_tree;
        reverse_needle_t _reverse_needle;
        uint32_t _error_count{};
    public:
        seed_prefix_extender(base_tree_t const & base_tree, needle_t needle, uint32_t error_count) noexcept :
            _base_tree{base_tree},
            _reverse_rcms{_base_tree.data().variants()},
            _reverse_tree{_reverse_rcms | libjst::make_volatile()},
            _reverse_needle{(needle_t &&) needle | std::views::reverse},
            _error_count{error_count}
        {}

        template <typename cargo_t, typename finder_t, typename callback_t>
        constexpr void operator()(cargo_t && seed_cargo,
                                  finder_t && seed_finder,
                                  callback_t && callback) const
        {
            // we need to map the position based on the local label!
            libjst::seek_position prefix_position = seed_cargo.position();
            if (std::ranges::empty(_reverse_needle)) {
                prefix_position.visit([&] <typename descriptor_t> (descriptor_t const &) {
                    if constexpr (!std::same_as<descriptor_t, libjst::breakpoint_end>)
                        prefix_position.initiate_alternate_node(prefix_position.get_variant_index());
                });
                callback(match_position{.tree_position = std::move(prefix_position),
                                        .label_offset = to_path_position(beginPosition(seed_finder), seed_cargo)},
                                        _error_count);
                return;
            }

            spm::restorable_myers_prefix_matcher extender{_reverse_needle, _error_count};
            auto reference_size = std::ranges::ssize(_reverse_tree.data().source());
            auto breakend_count = _reverse_tree.data().variants().size();

            // now reverse position!
            std::ptrdiff_t distance_to_end = std::ranges::ssize(seed_cargo.sequence()) - (beginPosition(seed_finder) - 1);
            std::ptrdiff_t global_start_offset = std::ranges::ssize(seed_cargo.path_sequence()) - distance_to_end;
            auto breakend_it = std::ranges::next(_base_tree.data().variants().begin(), prefix_position.get_variant_index());
            // TODO: what if breakend it high deletion breakend or low deletion breakend?
            // std::ptrdiff_t variant_position = libjst::position(*breakend_it);
            // bool prefix_starts_left_of_low = global_start_offset <= variant_position;

            match_position start{};

            start.label_offset = reference_size - global_start_offset;

            auto get_position = [] (auto const & breakend) {
                if (libjst::alt_kind(breakend) == libjst::alternate_sequence_kind::deletion)
                    return libjst::position(breakend);
                else
                    return libjst::high_breakend(breakend);
            };

            prefix_position.visit(libjst::multi_invocable{
                [&] (libjst::breakpoint_end) {
                    // log_debug("On reference path");
                    while (breakend_it != _base_tree.data().variants().begin() &&
                           global_start_offset < get_position(*breakend_it)) {
                        --breakend_it;
                    }
                    // log_debug("position forward low: ", get_position(*breakend_it)); // high position
                    ++breakend_it;
                    // log_debug("position forward high: ", libjst::position(*breakend_it)); // must be low position
                    assert(libjst::position(*breakend_it) >= global_start_offset);

                    std::ptrdiff_t low_breakend_index = std::ranges::distance(_base_tree.data().variants().begin(), breakend_it);
                    std::ptrdiff_t reverse_index_of_forward_low = breakend_count - low_breakend_index - 1;
                    start.tree_position.reset(reverse_index_of_forward_low, libjst::breakpoint_end::high);
                },
                [&] (libjst::alternate_path_descriptor const &) {
                    // log_debug("On alternate path");
                    if (global_start_offset >= libjst::position(*breakend_it)) {
                        // Then we set to the low boundary
                        // The low boundary is now an alternate node, and we have to enter this node in reverse direction.
                        std::ptrdiff_t reverse_index_of_forward_low = breakend_count - prefix_position.get_variant_index() - 1;
                        start.tree_position.initiate_alternate_node(reverse_index_of_forward_low);
                    } else {
                        while (breakend_it != _base_tree.data().variants().begin() &&
                               global_start_offset < get_position(*breakend_it)) {
                            --breakend_it;
                        }
                        // log_debug("position forward low: ", get_position(*breakend_it)); // high position
                        ++breakend_it;
                        // log_debug("position forward high: ", libjst::position(*breakend_it)); // must be low position
                        assert(libjst::position(*breakend_it) >= global_start_offset);

                        std::ptrdiff_t low_breakend_index = std::ranges::distance(_base_tree.data().variants().begin(), breakend_it);
                        std::ptrdiff_t reverse_index_of_forward_low = breakend_count - low_breakend_index - 1;
                        start.tree_position.reset(reverse_index_of_forward_low, libjst::breakpoint_end::high);
                    }

                }
            });
            // if (prefix_starts_left_of_low) {
            //     // we need to go back and find the first one that is not at this position.
            //     // TODO: Should it be the lowest ID?



            // } else {
            //     prefix_position.visit(libjst::multi_invocable{
            //         [&] (libjst::breakpoint_end site) {
            //             std::ptrdiff_t reverse_index_of_forward_low = breakend_count - prefix_position.get_variant_index();
            //             // The forward low represents the reverse high, so that the next forward breakend represents the
            //             // reverse low. Since the reverse index decreases with an increasing forward index we subtract one.
            //             // TODO: Check if going only to next. is sufficient or if we have to skip multiple if on the same position.
            //             start.tree_position.reset(reverse_index_of_forward_low - 1, libjst::breakpoint_end::high);
            //         },
            //         [&] (libjst::alternate_path_descriptor const &) {
            //             // The low boundary is now an alternate node, and we have to enter this node in reverse direction.
            //             std::ptrdiff_t reverse_index_of_forward_low = breakend_count - prefix_position.get_variant_index();
            //             start.tree_position.initiate_alternate_node(reverse_index_of_forward_low);
            //         }
            //     });
            //     // the global offset is inside of the variant or it must reset the node by one.
            // }

            // if global_start_offset >= variant position? => count - var_index ->
            // if (global_start_offset < variant_position) {
            //     log_debug("We can variant seek!");
            //     // we know the seed variant index is upper bound!
            //     // we get its index by subtracting
            //     start.tree_position.reset(breakend_count - seed_cargo.position().get_variant_index(),
            //                               libjst::breakpoint_end::low); // must site be the inverse?
            // } else {
            //     log_debug("We can use normal seek!");
            //     start.tree_position = seed_prefix_seek_position{seed_cargo.position(), breakend_count};
            // }

            // log_debug("reference_size: ", reference_size);
            // log_debug("global_start_offset: ", global_start_offset);
            // log_debug("low forward position: ", variant_position);
            // log_debug("prefix extender::start: ", start);
            // log_debug("prefix extender::start reverse: ", start.label_offset);
            // auto reverse_breakend_it = std::ranges::next(_reverse_rcms.variants().begin(), start.tree_position.get_variant_index());
            // log_debug("position reverse low: ", libjst::position(*reverse_breakend_it) - 1);
            // log_debug("position reverse high: ", libjst::position(*(reverse_breakend_it + 1)) - 1);

            // match_position start{.tree_position = seed_prefix_seek_position{seed_cargo.position(), breakend_count},
            //                      .label_offset = reference_size - (std::ranges::ssize(seed_cargo.path_sequence()) - distance_to_end)};
            // log_debug("prefix extender::cargo_sequence: ", std::ranges::ssize(seed_cargo.sequence()));
            // log_debug("prefix extender::cargo_path_sequence: ", std::ranges::ssize(seed_cargo.path_sequence()));
            // log_debug("prefix extender::finder_begin: ", beginPosition(seed_finder));
            // log_debug("prefix extender::finder_end: ", endPosition(seed_finder));
            // log_debug("prefix extender::global_begin: ", std::ranges::ssize(seed_cargo.path_sequence()) - distance_to_end);
            // log_debug("prefix extender::distance_to_end: ", distance_to_end);

            // log_debug("position reverse explicit: ", reference_size - variant_position);



            auto extend_tree = _reverse_tree | libjst::labelled()
                                             | libjst::coloured()
                                             | libjst::prune()
                                             | libjst::merge()
                                             | libjst::seek()
                                             | jstmap::extend_from(start, spm::window_size(extender));

            libjst::tree_traverser_base prefix_traverser{extend_tree};
            extension_state_manager manager{extender};
            prefix_traverser.subscribe(manager);
            // size_t counter{};
            for (auto cargo : prefix_traverser) {
                seed_prefix_node_cargo prefix_cargo{std::move(cargo), _reverse_tree};
                extender(prefix_cargo.sequence(), [&] (auto const & prefix_finder) {
                    // ++counter;
                    auto [best_position, best_score] = manager.top().second;
                    if (int32_t score = getScore(extender.capture()); score > best_score) {
                        manager.top().second =
                            std::pair{match_position{.tree_position = prefix_cargo.position(),
                                                     .label_offset = reference_size -
                                                        to_path_position(endPosition(prefix_finder), prefix_cargo)},
                                                    score};
                    }
                });
                if (cargo.is_leaf()) {
                    auto [best_position, best_score] = manager.top().second;
                    if (-best_score <= static_cast<decltype(best_score)>(_error_count))
                        callback(std::move(best_position), -best_score);
                }
            }
        }

    private:
        template <typename position_t, typename cargo_t>
        constexpr auto to_path_position(position_t local_position, cargo_t const & cargo) const noexcept {
            return std::ranges::ssize(cargo.path_sequence()) -
                   (std::ranges::ssize(cargo.sequence()) - local_position);
        }
    };
}  // namespace jstmap
