// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides pigeonhole filtration.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libspm/matcher/pigeonhole_matcher.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

namespace jstmap
{

    template <typename bucket_t>
    class pigeonhole_filter
    {
    private:

        bucket_t const & _bucket;
        double _error_rate{};

    public:

        explicit pigeonhole_filter(bucket_t const & bucket, double error_rate) noexcept :
            _bucket{bucket},
            _error_rate{error_rate}
        {}

        template <typename callback_t>
        constexpr void operator()(callback_t && callback) const {
            spm::pigeonhole_matcher filter{_bucket.needle_list, _error_rate};

            assert(spm::window_size(filter) > 0);
            // std::cout << "spm::window_size(filter) = " << spm::window_size(filter) << "\n";

            auto filter_tree = _bucket.base_tree | libjst::labelled()
                                                 | libjst::coloured()
                                                 | libjst::trim(spm::window_size(filter) - 1)
                                                 | libjst::prune()
                                                 | libjst::left_extend(spm::window_size(filter) - 1)
                                                 | libjst::merge()
                                                 | libjst::seek();

            libjst::tree_traverser_base traverser{filter_tree};
            for (auto it = traverser.begin(); it != traverser.end(); ++it) {
                auto seed_cargo = *it;
                // std::cout << "First sequence\n";
                filter(seed_cargo.sequence(), [&] (auto const & seed_finder) {
                    // std::cout << "Found seed\n";
                    callback(seed_cargo, seed_finder, filter.position());
                });
            }
        }
    };
}  // namespace jstmap
