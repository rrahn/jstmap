// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides search match aligner.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

#include <libspm/std/tag_invoke.hpp>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/match_position.hpp>
#include <jstmap/global/search_match.hpp>
#include <jstmap/search/type_alias.hpp>

namespace jstmap
{

    namespace detail {
        namespace _defer {
            template <auto & cpo, typename ...args_t>
            struct _cpo  {
                constexpr auto operator()() const
                    noexcept(std::is_nothrow_invocable_v<std::tag_t<cpo>, args_t...>)
                    -> std::invoke_result_t<std::tag_t<cpo>, args_t...> {
                    return std::invoke(cpo, std::declval<args_t>()...);
                }
            };

            template <auto & cpo, typename ...args_t>
            inline constexpr _cpo<cpo, args_t...> defer;
        } // namespace _defer

        using _defer::defer;

        template <auto & cpo>
        using instantiate_t = std::invoke_result_t<std::tag_t<cpo>>;

    } // namespace detail

    template <typename rcs_store_t, auto & ...tree_cpos>
    using composed_tree_t =
        std::remove_cvref_t<
            decltype((std::declval<rcs_store_t>() | ... | std::declval<detail::instantiate_t<tree_cpos>>()))
        >;

    class match_aligner {

        using ref_tree_type = composed_tree_t<rcs_store_t const &,
                                                libjst::make_volatile,
                                                libjst::labelled,
                                                // libjst::coloured,
                                                detail::defer<libjst::trim, size_t>,
                                                // libjst::prune,
                                                detail::defer<libjst::left_extend, size_t>,
                                                libjst::merge,
                                                libjst::seek
                                             >;

        record_sequence_t const & _query_sequence;
        ref_tree_type _reference_tree;
        double _error_rate{};

    public:

        explicit match_aligner(rcs_store_t const & rcs_store,
                               record_sequence_t const & query_sequence,
                               double error_rate = 0.03) :
            _query_sequence{query_sequence},
            _reference_tree{init(rcs_store)},
            _error_rate{error_rate}
        {}

        search_match operator()(match_position pos) const;

    private:

        template <typename match_score_t, typename mismatch_score_t,
                  typename gap_open_t, typename gap_extension_t>
        static auto get_alignment_config(match_score_t ms, mismatch_score_t mms, gap_open_t go, gap_extension_t ge) noexcept {
            return  seqan3::align_cfg::method_global{
                        seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                        seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                        seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                        seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}
                    }
                    | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{ms,mms}}
                    | seqan3::align_cfg::gap_cost_affine{go, ge};
        }

        ref_tree_type init(rcs_store_t const & rcs_store);
    };
}  // namespace jstmap
