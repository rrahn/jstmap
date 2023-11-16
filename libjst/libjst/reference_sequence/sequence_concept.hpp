// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for sequence.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename object_t>
    concept sequence = std::ranges::random_access_range<object_t> && std::ranges::sized_range<object_t>;

}  // namespace libjst
