// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_event_shared.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <vector>

#include <cereal/types/base_class.hpp>

#include <seqan3/alphabet/concept.hpp>

#include <libjst/detail/delta_event.hpp>
#include <libjst/utility/bit_vector_adaptor.hpp>

namespace libjst::detail
{

/*!\brief A delta event augmented with the information of how many sequences share this event.
 * \tparam alphabet_t The type of the alphabet of the underlying sequences; must model seqan3::semialphabet.
 *
 * \details
 *
 * In the libjst::journaled_sequence_tree sequences are represented through a reference sequence and a list of
 * shared delta events. They use an additional bitvector to store which sequence covers a particular delta event.
 * From these information the original sequences can be reconstructed or the journal sequence tree can be traversed
 * in a more efficient way, by only visiting all unique context specific subsequences.
 */
template <seqan3::semialphabet alphabet_t>
class delta_event_shared : public delta_event<alphabet_t>
{
private:
    //!\brief The type of the base delta event.
    using base_event_t = delta_event<alphabet_t>;
public:

    /*!\name Associated Types
     * \{
     */
    using coverage_type = libjst::utility::bit_vector_adaptor; //!< The coverage type
    using delta_event_type = base_event_t; //!< The original event type.
    using typename base_event_t::substitution_type;
    using typename base_event_t::insertion_type;
    using typename base_event_t::deletion_type;
    using typename base_event_t::alphabet_type;
    using typename base_event_t::segment_type;
    using typename base_event_t::size_type;
    using typename base_event_t::delta_variant_type;
    //!\}

private:
    coverage_type _coverage{}; //!< The coverage for this delta event.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    delta_event_shared() = default; //!< Default.
    delta_event_shared(delta_event_shared const &) = default; //!< Default.
    delta_event_shared(delta_event_shared &&) = default; //!< Default.
    delta_event_shared & operator=(delta_event_shared const &) = default; //!< Default.
    delta_event_shared & operator=(delta_event_shared &&) = default; //!< Default.
    ~delta_event_shared() = default; //!< Default.

    /*!\brief Constructs the shared delta event from a delta event and the coverage
     *
     * \param[in] event The delta event to represent.
     * \param[in] coverage The associated event coverage.
     */
    delta_event_shared(base_event_t && event, coverage_type coverage) :
        base_event_t{std::move(event)},
        _coverage{std::move(coverage)}
    {}

    /*!\brief Constructs the shared delta event from a position, event kind, and the associated coverage.
     *
     * \param[in] position The position of the delta event.
     * \param[in] event_kind The kind of the delta event.
     * \param[in] coverage The associated event coverage.
     */
    delta_event_shared(size_type const position, delta_variant_type event_kind, coverage_type coverage) :
        delta_event_shared{delta_event_type{position, std::move(event_kind)}, std::move(coverage)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the associated coverage.
    constexpr coverage_type & coverage() noexcept
    {
        return _coverage;
    }

    //!\overload
    constexpr coverage_type const & coverage() const noexcept
    {
        return _coverage;
    }
    //\}

    /*!\name Comparison
     * \{
     */
    using base_event_t::operator==;
    using base_event_t::operator<=>;

    //!\brief Compares for equality with other delta events.
    bool operator==(delta_event_shared const &) const = default;
    //!\}

    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves this shared delta event to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in, out] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(cereal::base_class<base_event_t>(this), _coverage);
    }

    /*!\brief Loads this shared delta event from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in, out] archive The archive to serialise this object from.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(cereal::base_class<base_event_t>(this), _coverage);
    }
    //!\}
};

/*!\name Type deduction guide
 * \{
 */
//!\brief Deduces the shared delta type from the passed libjst::detail::delta_event.
//!\relates libjst::detail::delta_event_shared
template <typename alphabet_t, typename coverage_t>
delta_event_shared(delta_event<alphabet_t>, coverage_t) -> delta_event_shared<alphabet_t>;
//!\}

/*!\brief Formatted output operator for the libjst::detail::delta_event_shared.
 * \relates libjst::detail::delta_event_shared
 *
 * \param[in, out] stream The output stream to write to.
 * \param[in] event The shared delta event to print.
 *
 * \returns A reference to the given output stream.
 */
template <typename char_t, typename char_traits_t, typename alphabet_t>
inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                              delta_event_shared<alphabet_t> const & event)
{
    using delta_event_t = typename delta_event_shared<alphabet_t>::delta_event_type;
    stream << static_cast<delta_event_t const &>(event)
           << " ~ <";
              std::ranges::for_each(event.coverage(), [&] (bool bit) mutable { stream << static_cast<int>(bit); });
    stream << '>';
    return stream;
}
}  // namespace libjst::detail
