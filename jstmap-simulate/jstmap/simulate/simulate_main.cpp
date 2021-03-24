// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map simulateer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

// ./jstmap simulate ../test/api/data/in.fasta ../test/api/output/tmp.jst

#include <seqan3/std/filesystem>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <seqan3/io/sequence_file/input.hpp>

#include <jstmap/simulate/simulate_main.hpp>
#include <jstmap/simulate/options.hpp>

#include <random> // uniform_int_distribution, random_device, mt19937
#include <set> // set
#include <stdexcept> // invalid_argument

#include <seqan3/core/debug_stream.hpp>

namespace jstmap
{

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::gapped<seqan3::dna5>; // instead of dna5
    using legal_sequence_alphabet_type = seqan3::dna5;
};
aligned_sequence_t load_sequence(std::filesystem::path const & sequence_file)
{
    aligned_sequence_t sequence;
    seqan3::sequence_file_input<my_traits> fin{sequence_file.c_str()};
    auto it = fin.begin();
    if (it == fin.end())
    {
        throw std::invalid_argument("Input file is empty.");
    }
    seqan3::debug_stream << seqan3::get<seqan3::field::seq>(*it) << '\n';
    sequence = seqan3::get<seqan3::field::seq>(*it);
    return sequence;
}
std::set<size_t> random_positions(size_t length, size_t n)
{
    static std::uniform_int_distribution<size_t> distr{0, length-1};
    static std::random_device engine;
    static std::mt19937 noise{engine()};
    std::set<size_t> positions;
    while (positions.size() < n) {
        positions.insert(distr(noise));
    }
    return positions;
}
seqan3::gapped<seqan3::dna5> random_char()
{
    static std::uniform_int_distribution<short> distr{0, 3};
    static std::random_device engine;
    static std::mt19937 noise{engine()};

    return static_cast<seqan3::gapped<seqan3::dna5> >(seqan3::dna4{}.assign_rank(distr(noise)));
    // return seqan3::gapped<seqan3::dna5>{}.assign_char(seqan3::dna5(seqan3::dna4{}.assign_rank(distr(noise)));
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
int simulate_main(seqan3::argument_parser & simulate_parser)
{
    simulate_options options{};

    simulate_parser.add_positional_option(options.input_file,
                                       "The input file.",
                                       seqan3::input_file_validator{{"fa", "fasta"}});
    simulate_parser.add_positional_option(options.output_file,
                                       "The output file.",
                                       seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                     {"jst"}});
    simulate_parser.add_option(options.error_rate,
                                       'e',
                                       "error-rate",
                                       "The relative error rate.");

    try
    {
        simulate_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    // Load the sequences.
    try
    {
        std::cout << "Loading sequences\n";
        auto sequence = load_sequence(options.input_file);
        aligned_sequence_t reference{sequence};
        // std::cout << "Actually not. Just creating some dummy sequences myself\n";
        // aligned_sequence_t reference(20);
        // for(size_t i = 0; i < reference.size(); ++i){
        //     reference[i] = static_cast<seqan3::dna5>(seqan3::dna4{}.assign_rank(i % 4));
        // }
        alignment_t alignment(reference, reference);

        std::set positions = random_positions(reference.size(), reference.size()*options.error_rate);
        for(auto it = positions.begin(); it != positions.end(); ++it)
        {
            std::cout << *it << " ";
        }
        std::cout << '\n';
        // // Choose n*error_rate*0.5 positions in alignment for SNPS
        // // Choose n*error_rate*0.25 positions in alignment for Insert
        // // Choose n*error_rate*0.25 positions in alignment for Delete
        // // iterator over both sequences
        // // iterator over three index vectors
        for(size_t i = 0; i < alignment.first.size(); ++i){
            std::cout << i%10 << " ";
        }
        std::cout << '\n';
        for(size_t i = 0; i < alignment.first.size(); ++i){
            seqan3::debug_stream << alignment.first[i] << " ";
        }
        std::cout << '\n';
        for(size_t i = 0; i < alignment.second.size(); ++i){
            seqan3::debug_stream << alignment.second[i] << " ";
        }


        std::cout << "\n\n";
        // if index of sequence iterator = index in SNP
            // exchange base in simulated by random different base
        // auto iterator_substitute = positions.begin();
        // for(size_t i = 0; i < alignment.first.size() && iterator_substitute != positions.end(); ++i){
        //     if (i == *iterator_substitute) {
        //         alignment.second[i] = (random_char(alignment.second[i]));
        //         ++iterator_substitute;
        //     }
        // }
        // if index of sequence iterator = index in Insert
            // add gap in reference
        // auto iterator_insert = positions.begin();
        // size_t j = 0;
        // for(size_t i = 0; i + j < alignment.first.size() && iterator_insert != positions.end(); ++i){
        //     if (i == *iterator_insert) {
        //         insert_gap(alignment.first,alignment.first.begin() + i + j);
        //         alignment.second.insert(alignment.second.begin() + i + j, random_char());
        //         ++iterator_insert; ++j;
        //         std::cout << i << " " << j << '\n';
        //     }
        // }
        // if index of sequence iterator = index in Delete
            // replace base with gap in reference
        auto iterator_delete = positions.begin();
        for(size_t i = 0; i < alignment.first.size() && iterator_delete != positions.end(); ++i){
            if (i == *iterator_delete) {
                alignment.second[i].assign_char('-');
                ++iterator_delete;
            }
        }


        for(size_t i = 0; i < alignment.first.size(); ++i){
            seqan3::debug_stream << alignment.first[i] << " ";
        }
        std::cout << '\n';
        for(size_t i = 0; i < alignment.second.size(); ++i){
            seqan3::debug_stream << alignment.second[i] << " ";
        }
        std::cout << '\n';
        // // auto tree = build_journaled_sequence_tree(std::move(sequences));
        // // serialise_jst(tree, options.output_file);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
