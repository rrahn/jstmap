#include <catch2/catch_test_macros.hpp>

#include <libjst/reference_sequence/sequence_breakpoint_simple.hpp>
#include <libjst/sequence_tree/breakpoint_sequence_label.hpp>

SCENARIO("Successfully creating a breakpoint sequence label", "[breakpoint_sequence_label]")
{
    GIVEN("a non-empty sequence and a simple breakpoint initialized with any positions")
    {
        std::string sequence{"is a sequence"};
        libjst::sequence_breakpoint_simple breakpoint{2, 5};

        WHEN("a breakpoint sequence label is instantiated with these parameters")
        {
            libjst::breakpoint_sequence_label label{sequence, breakpoint};

            THEN("the member function sequence returns the passed sequence")
            {
                REQUIRE(label.sequence() == sequence);
            }

            AND_THEN("the non-member function low_breakend returns the low_breakend of the passed breakpoint")
            {
                REQUIRE(libjst::low_breakend(label) == libjst::low_breakend(breakpoint));
                REQUIRE(libjst::low_breakend(std::as_const(label)) == libjst::low_breakend(breakpoint));
                REQUIRE(libjst::low_breakend(libjst::breakpoint_sequence_label{label}) ==
                        libjst::low_breakend(breakpoint));
            }

            AND_THEN("the non-member function high_breakend returns the high_breakend of the passed breakpoint")
            {
                REQUIRE(libjst::high_breakend(label) == libjst::high_breakend(breakpoint));
                REQUIRE(libjst::high_breakend(std::as_const(label)) == libjst::high_breakend(breakpoint));
                REQUIRE(libjst::high_breakend(libjst::breakpoint_sequence_label{label}) ==
                        libjst::high_breakend(breakpoint));
            }

            AND_THEN("the breakend_span is equal to the breakend_span of the passed breakpoint")
            {
                REQUIRE(libjst::breakend_span(label) == libjst::breakend_span(breakpoint));
                REQUIRE(libjst::breakend_span(std::as_const(label)) == libjst::breakend_span(breakpoint));
                REQUIRE(libjst::breakend_span(libjst::breakpoint_sequence_label{label}) ==
                        libjst::breakend_span(breakpoint));
            }
        }
    }
}
