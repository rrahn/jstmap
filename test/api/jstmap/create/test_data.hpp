// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap::test
{

using spm::operator""_dna4;

inline const raw_sequence_t reference = "aacctt"_dna4;

inline const std::vector<raw_sequence_t> sequences
{
    "aaaaaa"_dna4,
    "cccccc"_dna4,
    "tttttt"_dna4
};

}  // namespace jstmap::test
