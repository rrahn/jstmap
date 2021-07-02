// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map simulator.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */


#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/global/jstmap_type_alias.hpp>
#include <jstmap/create/serialise_jst.hpp>
#include <jstmap/simulate/load_reference.hpp>
#include <jstmap/simulate/options.hpp>
#include <jstmap/simulate/simulate_alignment.hpp>
#include <jstmap/simulate/simulate_main.hpp>

namespace jstmap
{

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
                               "The relative error rate.",
                               seqan3::option_spec::standard,
                               seqan3::arithmetic_range_validator{0,1});

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
        raw_sequence_t reference = load_reference(options.input_file);
        auto simulated = simulate_alignment(reference, options.error_rate);
        jst_t jst{std::move(reference)};
        jst.add(simulated);
        partitioned_jst_t partitioned_jst(std::addressof(jst), options.bin_count);
        std::cout << "Serialising JST\n";
        serialise(jst, partitioned_jst, options.output_file);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
