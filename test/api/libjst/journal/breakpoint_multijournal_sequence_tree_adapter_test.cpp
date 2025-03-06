// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------
#include <catch2/catch_test_macros.hpp>

#include <libjst/journal/breakpoint_multijournal.hpp>

SCENARIO("A breakpoint multijournal can be converted to a breakpoint_sequence_tree", "[breakpoint_multijournal_sequence_tree_adapter]")
{
    GIVEN("a breakpoint multijournal initialized with a non-empty std::string and 5 generated breakpoints")
    {
        std::string source{"AAAACCCCCGGGGGTTTTT"};
        libjst::breakpoint_multijournal journal{source};

        using breakpoint_t = libjst::sequence_breakpoint_t<std::string const &>;

        journal.record(breakpoint_t{1, 4}, std::string{});
        journal.record(breakpoint_t{3, 3}, std::string{"IIIIII"});
        journal.record(breakpoint_t{10, 11}, std::string{"J"});
        journal.record(breakpoint_t{13, 16}, std::string{});
        journal.record(breakpoint_t{13, 14}, std::string{"K"});


        WHEN("calling to_sequence_tree on this instance")
        {
            auto bst = libjst::to_sequence_tree(journal);
            auto node = bst.root();

            THEN("the root of the tree covers the first slice of the source until the low breakend of the first stored record")
            {
                REQUIRE(node->sequence() == source.substr(0,1));
                REQUIRE(libjst::low_breakend(*node) == 0ul);
                REQUIRE(libjst::high_breakend(*node) == 1ul);
            }

            AND_THEN("the sink of the tree returns a sentinel which is reached after calling 6 times next_ref() on the root")
            {
                for (size_t i = 0; i < 6; ++i)
                    node = node.next_ref();

                REQUIRE(node == bst.sink());
            }
        }
    }
}

