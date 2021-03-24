// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <string_view>
#include <string>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>

#include <jstmap/index/index_main.hpp> // Pulls in the index sub-command.
#include <jstmap/search/search_main.hpp> // Pulls in the search sub-command.
#include <jstmap/simulate/simulate_main.hpp> // Pulls in the search sub-command.

namespace jstmap
{

struct tool_names
{
    inline static const std::string base{"jstmap"};
    inline static const std::string index{"index"};
    inline static const std::string search{"search"};
    inline static const std::string simulate{"simulate"};

    static std::string subparser_name_for(std::string_view subcommand)
    {
        return std::string{base.data()} + "-" + std::string{subcommand.data()};
    }
};

} // namespace jstmap

int main(int const argc, char * const argv[])
{

    seqan3::argument_parser jstmap_parser{jstmap::tool_names::base, argc, argv, seqan3::update_notifications::off,
                                          {jstmap::tool_names::index, jstmap::tool_names::search, jstmap::tool_names::simulate}};

    jstmap_parser.info.description.push_back("The famous population mapper based on journaled string trees.");

    // TODO: What is best strategy to handle exceptions?
    // TODO: We need a logger instance: Should be allowed to set verbosity.
    try
    {
        jstmap_parser.parse();

        seqan3::argument_parser & selected_parser = jstmap_parser.get_sub_parser();

        if (selected_parser.info.app_name == jstmap::tool_names::subparser_name_for(jstmap::tool_names::index))
            return jstmap::index_main(selected_parser);
        else if (selected_parser.info.app_name == jstmap::tool_names::subparser_name_for(jstmap::tool_names::search))
            return jstmap::search_main(selected_parser);
        else if (selected_parser.info.app_name == jstmap::tool_names::subparser_name_for(jstmap::tool_names::simulate))
            return jstmap::simulate_main(selected_parser);
        else
            std::cerr << "Unknown subparser: " << selected_parser.info.app_name << '\n';
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        std::cerr << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    // Note: Arriving in this else branch means you did not handle all sub_parsers in the if branches above.

    return 0;

}
