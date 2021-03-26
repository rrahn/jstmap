// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------


#pragma once

#include <jstmap/simulate/global_types.hpp>
#include <map> // map

namespace jstmap
{

std::map<size_t, short> random_positions(size_t length, size_t n);
seqan3::gapped<seqan3::dna5> random_char();
seqan3::gapped<seqan3::dna5> random_char(seqan3::gapped<seqan3::dna5> old_char);
alignment_t simulate_alignment(sequence_t & reference, double error_rate);

}
