// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the simulation function, that generates a random alignment for a reference sequence.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <map> // map
#include <math.h> // ceil
#include <random> // uniform_int_distribution, random_device, mt19937

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/simulate/simulate_alignment.hpp>

namespace jstmap
{

template <std::uniform_random_bit_generator random_generator_t>
std::map<size_t, short> generate_random_positions(size_t length, size_t n, random_generator_t & generator)
{
    std::uniform_int_distribution<size_t> distr{0, length-1};
    std::map<size_t, short> positions;
    short error_type = 0;
    while (positions.size() < n)
    {
        const auto [it, success] = positions.insert({distr(generator), error_type});
        error_type += success;
        error_type &= 3; // mod 4
    }
    return positions;
}

template <std::uniform_random_bit_generator random_generator_t>
seqan3::gapped<spm::dna5> random_char(random_generator_t & generator)
{
    std::uniform_int_distribution<short> distr{0, 3};
    std::array tmp{'A', 'C', 'G', 'T'};
    return static_cast<seqan3::gapped<spm::dna5> >(
        seqan3::assign_char_to(tmp[distr(generator)], spm::dna5{}));
}

template <std::uniform_random_bit_generator random_generator_t>
seqan3::gapped<spm::dna5> random_char(seqan3::gapped<spm::dna5> old_char, random_generator_t & generator)
{
    std::uniform_int_distribution<short> distr{0, 3};
    std::array tmp{'A', 'C', 'G', 'T'};
    seqan3::gapped<spm::dna5> new_char;
    do
    {
        new_char = seqan3::assign_char_to(tmp[distr(generator)], spm::dna5{});
    } while (new_char == old_char);
    return new_char;
}

alignment_t simulate_alignment(raw_sequence_t & unaligned, double error_rate)
{
    assert(error_rate >= 0.0);
    assert(error_rate <= 1);

    aligned_sequence_t aligned{};
    seqan3::assign_unaligned(aligned, unaligned);
    alignment_t alignment{aligned, aligned};

    std::random_device engine;
    std::mt19937 noise{engine()};
    std::map positions = generate_random_positions(aligned.size(), ceil(aligned.size() * error_rate), noise);
    size_t j = 0;
    for(auto it = positions.begin(); it != positions.end(); ++it)
    {
        if (it->second == 3)
        {
            alignment.second[it->first + j].assign_char('-');
        }
        else if (it->second == 2)
        {
            seqan3::insert_gap(alignment.first, alignment.first.begin() + it->first + j);
            alignment.second.insert(alignment.second.begin() + it->first + j, random_char(noise));
            ++j;
        }
        else
        {
            alignment.second[it->first + j] = (random_char(alignment.second[it->first + j], noise));
        }
    }
    return alignment;
}

} // namespace jstmap
