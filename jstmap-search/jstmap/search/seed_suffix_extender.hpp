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

#include <stack>

#include <libspm/matcher/myers_prefix_matcher_restorable.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/partial_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/match_position.hpp>
#include <jstmap/search/seed_node_wrapper.hpp>
#include <jstmap/search/seed_extension_tree.hpp>

namespace jstmap
{
    template <typename base_tree_t, typename needle_t>
    class seed_suffix_extender
    {
        base_tree_t const & _base_tree;
        needle_t _needle;
        uint32_t _error_count{};
    public:
        seed_suffix_extender(base_tree_t const & base_tree, needle_t needle, uint32_t error_count) noexcept :
            _base_tree{base_tree},
            _needle{(needle_t &&) needle},
            _error_count{error_count}
        {}

        template <typename seed_cargo_t, typename finder_t, typename callback_t>
        constexpr void operator()(seed_cargo_t && seed_cargo,
                                  finder_t && seed_finder,
                                  callback_t && callback) const
        {
            if (std::ranges::empty(_needle)) {
                callback(match_position{.tree_position = seed_cargo.position(),
                                        .label_offset = endPosition(seed_finder)},
                                        _error_count);
                return;
            }

            spm::restorable_myers_prefix_matcher extender{_needle, _error_count};

            std::ptrdiff_t distance_to_end = std::ranges::ssize(seed_cargo.sequence()) - endPosition(seed_finder);
            match_position start{.tree_position = seed_cargo.position(),
                                 .label_offset = std::ranges::ssize(seed_cargo.path_sequence()) - distance_to_end};
            // log_info("Generate tree at: ", start);

            auto extend_tree = _base_tree | libjst::labelled()
                                          | libjst::coloured()
                                          | libjst::prune()
                                          | libjst::merge()
                                          | libjst::seek()
                                          | jstmap::extend_from(start, spm::window_size(extender));

            libjst::tree_traverser_base suffix_traverser{extend_tree};
            extension_state_manager manager{extender};
            suffix_traverser.subscribe(manager);
            // size_t counter{};
            for (auto cargo : suffix_traverser) {
                extender(cargo.sequence(), [&] (auto && suffix_finder) {
                    auto [best_position, best_score] = manager.top().second;
                    if (int32_t score = getScore(extender.capture()); score > best_score) {
                        manager.top().second = std::pair{match_position{.tree_position = cargo.position(),
                                                                        .label_offset = endPosition(suffix_finder)},
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
    };
}  // namespace jstmap
