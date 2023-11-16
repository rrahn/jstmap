// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides to_sequence_tree CPO.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/utility/tag_invoke.hpp>

namespace libjst
{

namespace _to_sequence_tree
{
inline constexpr struct _tag
{
    template <typename object_t>
        requires tag_invocable<_tag, object_t>
    constexpr auto operator()(object_t &&object) const noexcept(nothrow_tag_invocable<_tag, object_t>)
        -> tag_invoke_result_t<_tag, object_t>
    {
        return tag_invoke(_tag{}, (object_t &&)object);
    }
} to_sequence_tree{};

} // namespace _to_sequence_tree

using _to_sequence_tree::to_sequence_tree;

} // namespace libjst
