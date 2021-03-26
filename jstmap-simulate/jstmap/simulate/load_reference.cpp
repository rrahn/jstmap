// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map simulateer.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <seqan3/std/filesystem>
#include <seqan3/io/sequence_file/input.hpp>
#include <jstmap/simulate/load_reference.hpp>
#include <stdexcept> // invalid_argument

namespace jstmap
{

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    // using sequence_alphabet = seqan3::dna5; // instead of dna5
    using legal_sequence_alphabet_type = seqan3::dna5;
};
sequence_t load_reference(std::filesystem::path const & sequence_file)
{
    sequence_t sequence;
    seqan3::sequence_file_input<my_traits> fin{sequence_file.c_str()};
    auto it = fin.begin();
    if (it == fin.end())
    {
        throw std::invalid_argument("Input file is empty.");
    }
    sequence = seqan3::get<seqan3::field::seq>(*it);
    return sequence;
}

} // namespace jstmap
