// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan/sequence.h>

#include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap
{

using bin_sequence_t = std::views::all_t<raw_sequence_t const &>;
using bin_t = seqan::StringSet<bin_sequence_t>;

} // namespace jstmap