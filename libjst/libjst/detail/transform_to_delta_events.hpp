// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::transform_to_delta_events.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/views/zip.hpp>

#include <libjst/detail/delta_event.hpp>

namespace libjst::detail
{

/*!\brief Transforms a given alignment into a collection of delta events.
 * \tparam alphabet_t The type of the alphabet; must model seqan3::semialphabet.
 * \tparam alignment_t The type of the alignment.
 *
 * \returns A vector over libjst::detail::delta_event s.
 *
 * \details
 *
 * Iterates over the alignment and transforms all differences as a collection of libjst::detail::delta_event s.
 * The algorithm joins contiguous differences of the same kind (substitution, insertion, or deletion) into a single
 * delta event starting at the first detected position. The first sequence is assumed to be the reference sequence
 * and the second sequence is the target sequence. Thus, the positions of the delta events are based on the location
 * of the difference with respect to the reference sequence. In addition, a gap in the reference sequence is transformed
 * into an insertion event and a gap in the target sequence is transformed into a deletion event.
 */
template <seqan3::semialphabet alphabet_t, typename alignment_t>
constexpr auto transform_to_delta_events(alignment_t const & alignment)
{
    using delta_event_t = delta_event<alphabet_t>;
    using substitution_t = typename delta_event_t::substitution_type;
    using insertion_t = typename delta_event_t::insertion_type;
    using deletion_t = typename delta_event_t::deletion_type;

    std::vector<delta_event_t> result{};

    auto && [reference, target] = alignment;

    // Extracts the sequence from the zipped range with gapped characters.
    // It assumes that there is no gap character inside.
    auto extract_inserted_sequence = [] (auto && zipped_range)
    {
        std::vector<alphabet_t> sequence{};
        sequence.resize(std::ranges::distance(zipped_range));
        std::ranges::copy(zipped_range | std::views::transform([] (auto && gapped_pair)
        {
            assert(std::get<1>(gapped_pair).template is_alternative<alphabet_t>());

            return std::get<1>(gapped_pair).template convert_to<alphabet_t>();
        }), sequence.begin());

        return sequence;
    };

    auto zip = seqan3::views::zip(reference, target);
    auto it = zip.begin();
    size_t reference_position{0};
    while (it != zip.end())
    {
        auto [reference_value, target_value] = *it;
        if (reference_value == seqan3::gap{}) // insertion
        {
            auto next_it = std::ranges::find_if(it, zip.end(), [] (auto value_pair)
            {
                return value_pair.first != seqan3::gap{};
            });

            result.emplace_back(reference_position,
                                insertion_t{extract_inserted_sequence(std::ranges::subrange{it, next_it})});
            it = next_it;
        }
        else if (target_value == seqan3::gap{}) // deletion
        {
            auto next_it = std::ranges::find_if(it, zip.end(), [] (auto value_pair)
            {
                return value_pair.second != seqan3::gap{};
            });

            size_t const deletion_size = std::ranges::distance(it, next_it);
            result.emplace_back(reference_position, deletion_t{deletion_size});
            reference_position += deletion_size;
            it = next_it;
        }
        else if (reference_value != target_value) // substitution
        {
            auto next_it = std::ranges::find_if(it, zip.end(), [] (auto value_pair)
            {
                return value_pair.first == value_pair.second ||
                       value_pair.first == seqan3::gap{} ||
                       value_pair.second == seqan3::gap{};
            });

            auto sequence = extract_inserted_sequence(std::ranges::subrange{it, next_it});
            size_t const substitution_size = std::ranges::size(sequence);
            result.emplace_back(reference_position, substitution_t{std::move(sequence)});
            reference_position += substitution_size;
            it = next_it;
        }
        else
        {
            ++it;
            ++reference_position;
        }
    }

    return result;
}

}  // namespace libjst::detail
