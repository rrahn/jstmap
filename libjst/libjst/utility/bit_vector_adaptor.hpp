// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector_adaptor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <vector>

namespace libjst::utility
{

/*!\brief A bit vector adaptor for the std::vector<bool>.
 *
 * \details
 *
 * Wraps a std::vector<bool> and augments it with element wise bit-operations.
 */
class bit_vector_adaptor : public std::vector<bool>
{
private:
    //!\brief The base type.
    using base_t = std::vector<bool>;

public:

    using base_t::base_t;

    //!\brief Bitwise and.
    bit_vector_adaptor & operator&=(bit_vector_adaptor const & rhs) noexcept
    {
        assert(size() == rhs.size());

        std::transform(begin(), end(),
                       rhs.begin(),
                       begin(),
                       [] (bool const bit1, bool const bit2) { return bit1 & bit2; });
        return *this;
    }

    //!\brief Bitwise or.
    bit_vector_adaptor & operator|=(bit_vector_adaptor const & rhs) noexcept
    {
        assert(size() == rhs.size());

        std::transform(begin(), end(),
                       rhs.begin(),
                       begin(),
                       [] (bool const bit1, bool const bit2) { return bit1 | bit2; });
        return *this;
    }

    //!\brief Bitwise xor.
    bit_vector_adaptor & operator^=(bit_vector_adaptor const & rhs) noexcept
    {
        assert(size() == rhs.size());

        std::transform(begin(), end(),
                       rhs.begin(),
                       begin(),
                       [] (bool const bit1, bool const bit2) { return bit1 ^ bit2; });
        return *this;
    }

    //!\brief Bitwise not.
    bit_vector_adaptor operator~() const noexcept
    {
        bit_vector_adaptor tmp{*this};

        std::transform(begin(), end(), tmp.begin(), [] (bool const bit) { return !bit; });
        return tmp;
    }

    //!\brief Tests if all bits are 1.
    bool all() const noexcept
    {
        return std::ranges::all_of(*this, [] (bool const bit) { return bit == true; });
    }

    //!\brief Tests if any bit is 1.
    bool any() const noexcept
    {
        return std::ranges::any_of(*this, [] (bool const bit) { return bit == true; });
    }

    //!\brief Tests if all bits are 0.
    bool none() const noexcept
    {
        return std::ranges::none_of(*this, [] (bool const bit) { return bit == true; });
    }

private:
    //!\brief Cast to the base type.
    base_t const & as_base() const noexcept
    {
        return static_cast<base_t const &>(*this);
    }

    //!\overload
    base_t & as_base() noexcept
    {
        return static_cast<base_t &>(*this);
    }
};

/*!\name Bit operations
 * \{
 */
//!\brief Bitwise and.
//!\relates libjst::utility::bit_vector_adaptor
inline bit_vector_adaptor operator&(bit_vector_adaptor const & lhs, bit_vector_adaptor const & rhs)
{
    assert(lhs.size() == rhs.size());

    bit_vector_adaptor tmp{lhs};
    return tmp &= rhs;
}

//!\brief Bitwise or.
//!\relates libjst::utility::bit_vector_adaptor
inline bit_vector_adaptor operator|(bit_vector_adaptor const & lhs, bit_vector_adaptor const & rhs)
{
    assert(lhs.size() == rhs.size());

    bit_vector_adaptor tmp{lhs};
    return tmp |= rhs;
}

//!\brief Bitwise xor.
//!\relates libjst::utility::bit_vector_adaptor
inline bit_vector_adaptor operator^(bit_vector_adaptor const & lhs, bit_vector_adaptor const & rhs)
{
    assert(lhs.size() == rhs.size());

    bit_vector_adaptor tmp{lhs};
    return tmp ^= rhs;
}
//!\}

}  // namespace libjst::utility
