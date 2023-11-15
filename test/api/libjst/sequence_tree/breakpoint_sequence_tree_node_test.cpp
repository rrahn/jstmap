#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <iostream>
#include <string>

#include <libjst/journal/breakpoint_multijournal.hpp>
#include <libjst/sequence_tree/breakpoint_sequence_tree_node.hpp>

SCENARIO("Traversing a breakpoint multijournal using breakpoint_sequence_tree_node", "[breakpoint_sequence_tree_node]")
{
    GIVEN("a breakpoint_multijournal initialized with a non-empty std::string and a single record")
    {
        std::string journal_source{"AAAACCCCGGGGTTTT"};
        libjst::breakpoint_multijournal journal{journal_source};

        std::string alt_sequence{"NNNN"};
        auto it = journal.record(
            libjst::to_breakpoint(journal.source(), journal.source().begin() + 2, journal.source().begin() + 5),
            alt_sequence);

        WHEN("creating a breakpoint_sequence_tree_node covering the first slice of the journal source")
        {
            libjst::breakpoint_sequence_tree_node node{journal};

            THEN("the created node has a value")
            {
                REQUIRE(node.is_nil() == false);
            }
            AND_THEN("the node label is equal to the first slice of the journal source")
            {
                REQUIRE(std::ranges::equal(node.value().sequence(), std::string{"AA"}));
            }
            AND_THEN("the node label represents the breakpoint pointing to first record")
            {
                REQUIRE(libjst::low_breakend(node.value()) == 0ul);
                REQUIRE(libjst::high_breakend(node.value()) == libjst::low_breakend(*it));
            }

            WHEN("calling next_ref on node")
            {
                auto next_node = node.next_ref();

                THEN("the returned node has a value")
                {
                    REQUIRE(next_node.is_nil() == false);
                }
                AND_THEN("the node label is equal to the second slice of the journal source")
                {
                    REQUIRE(std::ranges::equal(next_node.value().sequence(), std::string{"AACCCCGGGGTTTT"}));
                }

                WHEN("calling next_alt on next_node")
                {
                    auto next_next_node = next_node.next_alt();

                    THEN("the returned node is nil")
                    {
                        REQUIRE(next_next_node.is_nil() == true);
                    }
                }
                AND_WHEN("calling next_ref on next_node")
                {
                    auto next_next_node = next_node.next_ref();

                    THEN("the returned node is nil")
                    {
                        REQUIRE(next_next_node.is_nil() == true);
                    }
                }
            }
            AND_WHEN("calling next_alt on node")
            {
                auto next_node = node.next_alt();

                THEN("the returned node has a value")
                {
                    REQUIRE(next_node.is_nil() == false);
                }
                AND_THEN("the node label is equal to the alternate sequence referenced by the inserted record")
                {
                    REQUIRE(std::ranges::equal(next_node.value().sequence(), alt_sequence));
                }

                WHEN("calling next_alt on next_node")
                {
                    auto next_next_node = next_node.next_alt();

                    THEN("the returned node is nil")
                    {
                        REQUIRE(next_next_node.is_nil() == true);
                    }
                }
                AND_WHEN("calling next_ref on next_node")
                {
                    auto next_next_node = next_node.next_ref();

                    THEN("the returned node is not nil")
                    {
                        REQUIRE(next_next_node.is_nil() == false);
                    }
                    AND_THEN("the label sequence is equal to the suffix of the journal source starting at 5")
                    {
                        REQUIRE(std::ranges::equal(next_next_node.value().sequence(), journal_source.substr(5)));
                    }
                }
            }
        }
    }
}

SCENARIO("Traversing a breakpoint multijournal with two insertion records at the same breakpoint",
         "[breakpoint_sequence_tree_node]")
{
    GIVEN("a breakpoint_multijournal initialized with a non-empty std::string and two insertion records at the same "
          "breakpoint")
    {
        std::string journal_source{"AAAACCCCGGGGTTTT"};
        libjst::breakpoint_multijournal journal{journal_source};

        std::string insertion1{"XX"};
        std::string insertion2{"YYY"};
        auto breakpoint =
            libjst::to_breakpoint(journal.source(), journal.source().begin() + 2, journal.source().begin() + 2);
        journal.record(breakpoint, insertion1);
        journal.record(breakpoint, insertion2);

        libjst::breakpoint_sequence_tree_node root{journal};
        WHEN("current node is root")
        {
            THEN("root represents first slice and breakends")
            {
                REQUIRE(std::ranges::equal(root.value().sequence(), std::string{"AA"}));
                REQUIRE(libjst::low_breakend(root.value()) == 0ul);
                REQUIRE(libjst::high_breakend(root.value()) == 2ul);
            }
            WHEN("calling next_ref() on root") {
                THEN("child returns node representing empty slice at [2,2)") {
                    auto child = root.next_ref();
                    REQUIRE(std::string{child.value().sequence().begin(), child.value().sequence().end()} == std::string{});
                    REQUIRE(libjst::low_breakend(child.value()) == 2ul);
                    REQUIRE(libjst::high_breakend(child.value()) == 2ul);
                }
            }
            AND_WHEN("calling next_alt() on root")
            {
                auto child = root.next_alt();
                THEN("child node represents insertion2")
                {
                    REQUIRE(std::ranges::equal(child.value().sequence(), insertion2));
                    REQUIRE(libjst::low_breakend(child.value()) == 2ul);
                    REQUIRE(libjst::high_breakend(child.value()) == 2ul);
                }
                AND_THEN("calling next_alt() on child returns nil")
                {
                    REQUIRE(child.next_alt().is_nil() == true);
                }
                AND_THEN("calling next_ref() on child returns empty source slice")
                {
                    auto grandchild = child.next_ref();
                    REQUIRE(std::ranges::equal(grandchild.value().sequence(), std::string{""}));
                    REQUIRE(libjst::low_breakend(grandchild.value()) == 2ul);
                    REQUIRE(libjst::high_breakend(grandchild.value()) == 2ul);
                }
            }
        }
        WHEN("current node is reference node between two insertions")
        {
            auto node = root.next_alt().next_ref();
            THEN("calling next_alt() returns node representing insertion1")
            {
                auto child = node.next_alt();
                REQUIRE(std::ranges::equal(child.value().sequence(), insertion1));
                REQUIRE(libjst::low_breakend(child.value()) == 2ul);
                REQUIRE(libjst::high_breakend(child.value()) == 2ul);
                WHEN("calling next_alt() on child") {
                    THEN("grandchild is nil") {
                        REQUIRE(child.next_alt().is_nil() == true);
                    }
                    AND_THEN("calling next_ref() on child returns source slice starting at 2") {
                        auto grandchild = child.next_ref();
                        REQUIRE(std::ranges::equal(grandchild.value().sequence(), journal_source.substr(2)));
                        REQUIRE(libjst::low_breakend(grandchild.value()) == 2ul);
                        REQUIRE(libjst::high_breakend(grandchild.value()) == 16ul);
                    }
                }
            }
            AND_THEN("calling next_ref() on child returns source slice starting at 2")
            {
                auto child = node.next_ref();
                REQUIRE(std::ranges::equal(child.value().sequence(), journal_source.substr(2)));
                REQUIRE(libjst::low_breakend(child.value()) == 2ul);
                REQUIRE(libjst::high_breakend(child.value()) == 16ul);
            }
        }
    }
}
