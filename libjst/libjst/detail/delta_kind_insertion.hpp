// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_kind_insertion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <cereal/types/base_class.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/range/sequence.hpp>

#include <libjst/detail/delta_kind_base.hpp>

namespace libjst::detail
{
/*!\brief A delta event representing an insertion.
 *
 * \tparam alphabet_t The alphabet type used to store the inserted sequence; must model seqan3::semialphabet.
 */
template <seqan3::semialphabet alphabet_t>
class delta_kind_insertion : public delta_kind_base<std::vector<alphabet_t>>
{
private:
    using base_t = delta_kind_base<std::vector<alphabet_t>>; //!< The base type.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr delta_kind_insertion() = default; //!< Default.
    constexpr delta_kind_insertion(delta_kind_insertion const &) = default; //!< Default.
    constexpr delta_kind_insertion(delta_kind_insertion &&) = default; //!< Default.
    constexpr delta_kind_insertion & operator=(delta_kind_insertion const &) = default; //!< Default.
    constexpr delta_kind_insertion & operator=(delta_kind_insertion &&) = default; //!< Default.
    ~delta_kind_insertion() = default; //!< Default.

    /*!\brief Initialises an insertion with the given insertion sequence.
     *
     * \tparam sequence_t The sequence type; must model seqan3::sequence.
     *
     * \param[in] sequence The inserted sequence.
     */
    template <seqan3::sequence sequence_t>
    explicit constexpr delta_kind_insertion(sequence_t && sequence) : base_t{std::forward<sequence_t>(sequence)}
    {}
    //!\}

    //!\brief Compare against other insertions for equality.
    bool operator==(delta_kind_insertion const &) const = default;

    /*!\name Serialisation
     * \{
     */
    //!\copydoc libjst::detail::delta_kind_base::save
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(cereal::base_class<base_t>(this));
    }

    //!\copydoc libjst::detail::delta_kind_base::load
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(cereal::base_class<base_t>(this));
    }
    //!\}
};

/*!\name Type deduction guide
 * \{
 */
/*!\brief Deduces the alphabet type from the given sequence.
 * \relates libjst::detail::delta_kind_insertion
 */
template <seqan3::sequence sequence_t>
delta_kind_insertion(sequence_t) -> delta_kind_insertion<std::ranges::range_value_t<sequence_t>>;
//!\}
}  // namespace libjst::detail
