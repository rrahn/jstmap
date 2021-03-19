// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_kind_deletion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/detail/delta_kind_base.hpp>

namespace libjst::detail
{
//!\brief A delta event representing a deletion.
class delta_kind_deletion : public delta_kind_base<size_t>
{
public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr delta_kind_deletion() = default; //!< Default.
    constexpr delta_kind_deletion(delta_kind_deletion const &) = default; //!< Default.
    constexpr delta_kind_deletion(delta_kind_deletion &&) = default; //!< Default.
    constexpr delta_kind_deletion & operator=(delta_kind_deletion const &) = default; //!< Default.
    constexpr delta_kind_deletion & operator=(delta_kind_deletion &&) = default; //!< Default.
    ~delta_kind_deletion() = default; //!< Default.

    /*!\brief Initialises a deletion with the given deletion size.
     *
     * \param[in] deletion_size The size of the deletion.
     */
    explicit constexpr delta_kind_deletion(size_t const deletion_size) noexcept : delta_kind_base<size_t>{deletion_size}
    {}
    //!\}

    //!\brief Compare against other deletions for equality.
    bool operator==(delta_kind_deletion const &) const = default;
};

}  // namespace libjst::detail
