// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_kind_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <concepts>

#include <seqan3/range/concept.hpp>

namespace libjst::detail
{
//!\brief A delta event representing a deletion.
template <std::semiregular value_t>
class delta_kind_base
{
private:
    //!\brief The stored value.
    value_t _value{};
public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr delta_kind_base() = default; //!< Default.
    constexpr delta_kind_base(delta_kind_base const &) = default; //!< Default.
    constexpr delta_kind_base(delta_kind_base &&) = default; //!< Default.
    constexpr delta_kind_base & operator=(delta_kind_base const &) = default; //!< Default.
    constexpr delta_kind_base & operator=(delta_kind_base &&) = default; //!< Default.
    ~delta_kind_base() = default; //!< Default.

    /*!\brief Initialises the delta event kind with the given value.
     *
     * \param[in] value The value to store.
     */
    explicit constexpr delta_kind_base(value_t value) noexcept : _value{std::move(value)}
    {}

    /*!\brief Initialises the delta event kind with the given value if it is a seqan3::sequence.
     *
     * \param[in] sequence The sequence to store.
     */
    template <seqan3::sequence sequence_t>
        requires seqan3::sequence<value_t>
    explicit constexpr delta_kind_base(sequence_t && sequence) noexcept
    {
        if constexpr (std::assignable_from<value_t &, sequence_t &&>)
        {
            _value = std::move(sequence);
        }
        else
        {
            _value.reserve(std::ranges::distance(sequence));
            std::ranges::move(std::forward<sequence_t>(sequence), std::back_inserter(_value));
        }
    }
    //!\}

    //!\brief Returns the stored value.
    value_t const & value() const noexcept
    {
        return _value;
    }

    //!\brief Compare against other delta kinds for equality.
    bool operator==(delta_kind_base const &) const = default;

    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves the delta kind to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_value);
    }

    /*!\brief Loads the delta kind from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in] archive The archive to serialise this object from.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(_value);
    }
    //!\}
};

}  // namespace libjst::detail
