// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides matching strategy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/matcher/horspool_matcher.hpp>
#include <libjst/traversal/state_oblivious_traverser.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>

#include <jstmap/search/type_alias.hpp>

namespace jstmap
{

    class matching_operation {
    private:

    public:

        matching_operation() = default;

        template <typename callback_t>
        constexpr void operator()(haystack_type && haystack, search_queries_type const & bucket, callback_t && callback) const {

            std::ranges::for_each(bucket, [&] (search_query const & query) {
                jst::contrib::horspool_matcher matcher{query.value().sequence()};

                if (jst::contrib::window_size(matcher) == 0)
                    return;

                auto search_tree = haystack | libjst::labelled()
                                            | libjst::coloured()
                                            | libjst::trim(jst::contrib::window_size(matcher) - 1)
                                            | libjst::prune_unsupported()
                                            | libjst::left_extend(jst::contrib::window_size(matcher) - 1)
                                            | libjst::merge() // make big nodes
                                            | libjst::seek();

                libjst::tree_traverser_base oblivious_path{search_tree};
                for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
                    auto && cargo = *it;
                    matcher(cargo.sequence(), [&] ([[maybe_unused]] auto && label_finder) {
                        callback(query, match_position{.tree_position{cargo.position()},
                                                       .label_offset{std::ranges::ssize(cargo.sequence()) - seqan2::endPosition(label_finder)}});
                    });
                }
            });
        }
    };

}  // namespace jstmap
