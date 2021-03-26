// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map indexer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

// #include <seqan3/std/filesystem>

// #include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/index/serialise_jst.hpp>
#include <jstmap/index/load_sequence.hpp>
#include <jstmap/index/index_main.hpp>
#include <jstmap/index/journaled_sequence_tree_builder.hpp>
#include <jstmap/index/options.hpp>

namespace jstmap
{

int index_main(seqan3::argument_parser & index_parser)
{
    index_options options{};

    index_parser.add_positional_option(options.input_file,
                                       "The input file.",
                                       seqan3::input_file_validator{{"fa", "fasta"}});
    index_parser.add_positional_option(options.output_file,
                                       "The output file.",
                                       seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                     {"jst"}});

    try
    {
        index_parser.parse();
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
        auto sequences = load_sequences(options.input_file);

        auto tree = build_journaled_sequence_tree(std::move(sequences));
        serialise_jst(tree, options.output_file);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
