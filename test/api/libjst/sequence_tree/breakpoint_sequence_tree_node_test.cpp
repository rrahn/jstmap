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

SCENARIO("Finding the first non-overlapping breakpoint after a deletion", "[breakpoint_sequence_tree_node]")
{
    GIVEN("a breakpoint multijournal initialized with a non-empty source sequence")
    {
        std::string journal_source{"AAAACCCCGGGGTTTT"};
        libjst::breakpoint_multijournal journal{journal_source};

        using breakpoint_t = libjst::sequence_breakpoint_t<std::string const &>;

        WHEN ("two non-overlapping deletions are recorded")
        {
            breakpoint_t breakpoint1{2, 4};
            breakpoint_t breakpoint2{6, 10};
            journal.record(breakpoint1, std::string{});
            journal.record(breakpoint2, std::string{});

            THEN ("the first alt node represents breakpoint1")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt1 = node.next_alt();
                REQUIRE(alt1.value().sequence().empty());
                REQUIRE(libjst::low_breakend(alt1.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt1.value()) == 4ul);

                AND_THEN("the next ref node covers the journal source from 4 to 6")
                {
                    auto alt1_ref = alt1.next_ref();
                    REQUIRE(alt1_ref.value().sequence() == journal_source.substr(4, 2));
                    REQUIRE(libjst::low_breakend(alt1_ref.value()) == 4ul);
                    REQUIRE(libjst::high_breakend(alt1_ref.value()) == 6ul);

                    AND_THEN("the next alt node covers the second breakpoint")
                    {
                        auto alt1_ref_alt2 = alt1_ref.next_alt();
                        REQUIRE(alt1_ref_alt2.value().sequence().empty());
                        REQUIRE(libjst::low_breakend(alt1_ref_alt2.value()) == 6ul);
                        REQUIRE(libjst::high_breakend(alt1_ref_alt2.value()) == 10ul);
                    }
                }
            }
        }
        WHEN ("recording two overlapping deletions")
        {
            breakpoint_t breakpoint1{2, 8};
            breakpoint_t breakpoint2{6, 10};
            journal.record(breakpoint1, std::string{});
            journal.record(breakpoint2, std::string{});

            THEN ("there is no path including both breakpoints")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt1 = node.next_alt();
                REQUIRE(alt1.value().sequence().empty());
                REQUIRE(libjst::low_breakend(alt1.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt1.value()) == 8ul);

                auto alt1_ref = alt1.next_ref();
                REQUIRE(alt1_ref.value().sequence() == journal_source.substr(8));
                REQUIRE(libjst::low_breakend(alt1_ref.value()) == 8ul);
                REQUIRE(libjst::high_breakend(alt1_ref.value()) == 16ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 4));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 6ul);

                auto alt2 = node.next_alt();
                REQUIRE(alt2.value().sequence().empty());
                REQUIRE(libjst::low_breakend(alt2.value()) == 6ul);
                REQUIRE(libjst::high_breakend(alt2.value()) == 10ul);

                auto alt2_ref = alt2.next_ref();
                REQUIRE(alt2_ref.value().sequence() == journal_source.substr(10));
                REQUIRE(libjst::low_breakend(alt2_ref.value()) == 10ul);
                REQUIRE(libjst::high_breakend(alt2_ref.value()) == 16ul);
            }
        }
        WHEN ("recording multiple records with overlap at a breakend followed by another breakpoint")
        {
            breakpoint_t breakpoint1{2, 2};
            breakpoint_t breakpoint2{2, 3};
            breakpoint_t breakpoint3{2, 3};
            breakpoint_t breakpoint4{2, 8};
            breakpoint_t breakpoint5{4, 6};
            breakpoint_t breakpoint6{6, 10};
            breakpoint_t breakpoint7{7, 7};
            breakpoint_t breakpoint8{8, 9};

            std::string alt1{"ZZZ"};
            std::string alt2{"X"};
            std::string alt3{"Y"};
            std::string alt4{};
            std::string alt5{};
            std::string alt6{};
            std::string alt7{"I"};
            std::string alt8{"J"};

            journal.record(breakpoint1, alt1);
            journal.record(breakpoint2, alt2);
            journal.record(breakpoint3, alt3);
            journal.record(breakpoint4, alt4);
            journal.record(breakpoint5, alt5);
            journal.record(breakpoint6, alt6);
            journal.record(breakpoint7, alt7);
            journal.record(breakpoint8, alt8);

            THEN("the reference path consists of 9 nodes")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                REQUIRE(node.is_nil() == false);
                REQUIRE(node.value().sequence() == journal_source.substr(0, 2));
                REQUIRE(libjst::low_breakend(node.value()) == 0ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 2));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 4ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(4, 2));
                REQUIRE(libjst::low_breakend(node.value()) == 4ul);
                REQUIRE(libjst::high_breakend(node.value()) == 6ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(6, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 6ul);
                REQUIRE(libjst::high_breakend(node.value()) == 7ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(7, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 7ul);
                REQUIRE(libjst::high_breakend(node.value()) == 8ul);

                node = node.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(8));
                REQUIRE(libjst::low_breakend(node.value()) == 8ul);
                REQUIRE(libjst::high_breakend(node.value()) == 16ul);

                node = node.next_ref();
                REQUIRE(node.is_nil() == true);
                REQUIRE(node == libjst::breakpoint_sequence_tree_sentinel{});
                REQUIRE(libjst::breakpoint_sequence_tree_sentinel{} == node);
            }
            AND_THEN("There is a path covering breakpoint1, breakpoint2, breakpoint5, and breakpoint6")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt_node1 = node.next_alt();
                REQUIRE(alt_node1.value().sequence() == alt1);
                REQUIRE(libjst::low_breakend(alt_node1.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt_node1.value()) == 2ul);

                node = alt_node1.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                auto alt_node2 = node.next_alt();
                REQUIRE(alt_node2.value().sequence() == alt2);
                REQUIRE(libjst::low_breakend(alt_node2.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt_node2.value()) == 3ul);

                node = alt_node2.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(3, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 3ul);
                REQUIRE(libjst::high_breakend(node.value()) == 4ul);

                auto alt_node5 = node.next_alt();
                REQUIRE(alt_node5.value().sequence() == alt5);
                REQUIRE(libjst::low_breakend(alt_node5.value()) == 4ul);
                REQUIRE(libjst::high_breakend(alt_node5.value()) == 6ul);

                node = alt_node5.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(6, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 6ul);
                REQUIRE(libjst::high_breakend(node.value()) == 6ul);

                auto alt_node6 = node.next_alt();
                REQUIRE(alt_node6.value().sequence() == alt6);
                REQUIRE(libjst::low_breakend(alt_node6.value()) == 6ul);
                REQUIRE(libjst::high_breakend(alt_node6.value()) == 10ul);

                node = alt_node6.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(10));
                REQUIRE(libjst::low_breakend(node.value()) == 10ul);
                REQUIRE(libjst::high_breakend(node.value()) == 16ul);
            }
            AND_THEN("there is a path covering breakpoint1, breakpoint3, breakpoint5, and breakpoint6")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt_node1 = node.next_alt();

                node = alt_node1.next_ref().next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                auto alt_node3 = node.next_alt();
                REQUIRE(alt_node3.value().sequence() == alt3);
                REQUIRE(libjst::low_breakend(alt_node3.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt_node3.value()) == 3ul);

                node = alt_node3.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(3, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 3ul);
                REQUIRE(libjst::high_breakend(node.value()) == 4ul);

                auto alt_node5 = node.next_alt();
                REQUIRE(alt_node5.value().sequence() == alt5);
                REQUIRE(libjst::low_breakend(alt_node5.value()) == 4ul);
                REQUIRE(libjst::high_breakend(alt_node5.value()) == 6ul);
            }
            AND_THEN("there is a path covering breakpoint1, breakpoint4, and breakpoint8")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt_node1 = node.next_alt();
                node = alt_node1.next_ref().next_ref().next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 2ul);

                auto alt_node4 = node.next_alt();
                REQUIRE(alt_node4.value().sequence() == alt4);
                REQUIRE(libjst::low_breakend(alt_node4.value()) == 2ul);
                REQUIRE(libjst::high_breakend(alt_node4.value()) == 8ul);

                node = alt_node4.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(8, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 8ul);
                REQUIRE(libjst::high_breakend(node.value()) == 8ul);

                auto alt_node8 = node.next_alt();
                REQUIRE(alt_node8.value().sequence() == alt8);
                REQUIRE(libjst::low_breakend(alt_node8.value()) == 8ul);
                REQUIRE(libjst::high_breakend(alt_node8.value()) == 9ul);

                node = alt_node8.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(9));
                REQUIRE(libjst::low_breakend(node.value()) == 9ul);
                REQUIRE(libjst::high_breakend(node.value()) == 16ul);
            }
            AND_THEN("there is a path covering breakpoint1 and breakpoint5")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                auto alt_node1 = node.next_alt();
                node = alt_node1.next_ref().next_ref().next_ref().next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 2));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 4ul);

                auto alt_node5 = node.next_alt();
                REQUIRE(alt_node5.value().sequence() == alt5);
                REQUIRE(libjst::low_breakend(alt_node5.value()) == 4ul);
                REQUIRE(libjst::high_breakend(alt_node5.value()) == 6ul);

                node = alt_node5.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(6, 0));
                REQUIRE(libjst::low_breakend(node.value()) == 6ul);
                REQUIRE(libjst::high_breakend(node.value()) == 6ul);
            }
            AND_THEN("there is a path covering breakpoint5, breakpoint7, and breakpoint8")
            {
                libjst::breakpoint_sequence_tree_node node{journal};
                node = node.next_ref().next_ref().next_ref().next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(2, 2));
                REQUIRE(libjst::low_breakend(node.value()) == 2ul);
                REQUIRE(libjst::high_breakend(node.value()) == 4ul);

                auto alt_node5 = node.next_alt();
                REQUIRE(alt_node5.value().sequence() == alt5);
                REQUIRE(libjst::low_breakend(alt_node5.value()) == 4ul);
                REQUIRE(libjst::high_breakend(alt_node5.value()) == 6ul);

                node = alt_node5.next_ref().next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(6, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 6ul);
                REQUIRE(libjst::high_breakend(node.value()) == 7ul);

                auto alt_node7 = node.next_alt();
                REQUIRE(alt_node7.value().sequence() == alt7);
                REQUIRE(libjst::low_breakend(alt_node7.value()) == 7ul);
                REQUIRE(libjst::high_breakend(alt_node7.value()) == 7ul);

                node = alt_node7.next_ref();
                REQUIRE(node.value().sequence() == journal_source.substr(7, 1));
                REQUIRE(libjst::low_breakend(node.value()) == 7ul);
                REQUIRE(libjst::high_breakend(node.value()) == 8ul);

                auto alt_node8 = node.next_alt();
                REQUIRE(alt_node8.value().sequence() == alt8);
                REQUIRE(libjst::low_breakend(alt_node8.value()) == 8ul);
                REQUIRE(libjst::high_breakend(alt_node8.value()) == 9ul);
            }
        }
    }
}
