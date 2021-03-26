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

#include <jstmap/simulate/simulate_alignment.hpp>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <random> // uniform_int_distribution, random_device, mt19937
#include <math.h> // ceil

#include <seqan3/core/debug_stream.hpp>

namespace jstmap
{

std::map<size_t, short> random_positions(size_t length, size_t n)
{
    static std::uniform_int_distribution<size_t> distr{0, length-1};
    static std::random_device engine;
    static std::mt19937 noise{engine()};
    std::map<size_t, short> positions;
    short error_type = 0;
    while (positions.size() < n) {
        const auto [it, success] = positions.insert({distr(noise), error_type});
        error_type += success;
        error_type &= 3; // mod 4
    }
    return positions;
}
seqan3::gapped<seqan3::dna5> random_char()
{
    static std::uniform_int_distribution<short> distr{0, 3};
    static std::random_device engine;
    static std::mt19937 noise{engine()};

    return static_cast<seqan3::gapped<seqan3::dna5> >(seqan3::dna4{}.assign_rank(distr(noise)));
}
seqan3::gapped<seqan3::dna5> random_char(seqan3::gapped<seqan3::dna5> old_char)
{
    static std::uniform_int_distribution<short> distr{0, 3};
    static std::random_device engine;
    static std::mt19937 noise{engine()};
    seqan3::gapped<seqan3::dna5> new_char;
    do {
        new_char = static_cast<seqan3::dna5>(seqan3::dna4{}.assign_rank(distr(noise)));
    } while (new_char == old_char);
    return new_char;
}
alignment_t simulate_alignment(aligned_sequence_t & reference, double error_rate)
{
    seqan3::arithmetic_range_validator error_validator{0,1};
    error_validator(error_rate);
    alignment_t alignment(reference, reference);
    for(size_t i = 0; i < alignment.first.size(); ++i){
        seqan3::debug_stream << alignment.first[i];
    }
    seqan3::debug_stream << '\n';
    for(size_t i = 0; i < alignment.second.size(); ++i){
        seqan3::debug_stream << alignment.second[i];
    }
    seqan3::debug_stream << '\n' << reference.size() << '\n' << ceil(reference.size()*error_rate) << '\n';
    std::map positions = random_positions(reference.size(), ceil(reference.size()*error_rate));
    size_t j = 0;
    for(auto it = positions.begin(); it != positions.end(); ++it){
        if (it->second == 3) {
            alignment.second[it->first + j].assign_char('-');
            seqan3::debug_stream << __LINE__ << '\n';
        } else if (it->second == 2) {
            seqan3::insert_gap(alignment.first, alignment.first.begin() + it->first + j);
            alignment.second.insert(alignment.second.begin() + it->first + j, random_char());
            ++j;
            seqan3::debug_stream << __LINE__ << '\n';
        } else {
            alignment.second[it->first + j] = (random_char(alignment.second[it->first + j]));
            seqan3::debug_stream << __LINE__ << '\n';
        }
    }
    for(size_t i = 0; i < alignment.first.size(); ++i){
        seqan3::debug_stream << alignment.first[i];
    }
    seqan3::debug_stream << '\n';
    for(size_t i = 0; i < alignment.second.size(); ++i){
        seqan3::debug_stream << alignment.second[i];
    }
    seqan3::debug_stream << '\n';
    return alignment;
}

} // namespace jstmap
