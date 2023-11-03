#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <string>
#include <ranges>

#include <libjst/journal/basic_multisequence_journal.hpp>

namespace test
{
    template <std::ranges::forward_range t>
    struct debug_sequence
    {
        t value;

        template <typename char_t, typename char_traits_t>
        friend std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream, debug_sequence const & seq)
        {
            for (auto && c : seq.value)
                stream << c;
            return stream;
        }
    };
}

SCENARIO("Recording insertions", "[basic_multisequence_journal][insertion]") {
    using namespace std::string_literals;

    GIVEN("A basic_multisequence_journal initialized with a std::string source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::basic_multisequence_journal journal{source};

        auto i = GENERATE(Catch::Generators::values({0, 8, 16}));
        WHEN("Inserting a sequence at different positions") {
            std::string alt_sequence{"xxx"};
            auto breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + i);
            auto it = journal.record(breakpoint, alt_sequence);
            THEN ("it points to the recorded insertion") {
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequence));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == i);
            } AND_THEN ("the journal contains one element") {
                REQUIRE(journal.size() == 1);
            }
        }

        WHEN("Inserting three sequences at the same position with different lengths") {
            std::vector alt_sequences{"x"s, "xx"s, "xxx"s};
            auto breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + i);
            journal.record(breakpoint, alt_sequences[0]);
            journal.record(breakpoint, alt_sequences[1]);
            journal.record(breakpoint, alt_sequences[2]);
            THEN ("the journal contains three elements") {
                REQUIRE(journal.size() == 3);
            } AND_THEN ("the elements are ordered by insertion length in descending order") {
                auto it = journal.begin();
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequences[2]));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == i);
                REQUIRE(std::ranges::equal((++it)->sequence(), alt_sequences[1]));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == i);
                REQUIRE(std::ranges::equal((++it)->sequence(), alt_sequences[0]));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == i);
            }
        }
    }
}

SCENARIO("Recording deletions", "[basic_multisequence_journal][deletion]") {
    using namespace std::string_literals;

    GIVEN("A basic_multisequence_journal initialized with a std::string source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::basic_multisequence_journal journal{source};

        WHEN("Recording a deletion at different breakpoints") {
            auto i = GENERATE(Catch::Generators::values({0, 8, 16}));
            auto j = GENERATE(Catch::Generators::values({1, 8, 16}));

            auto breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + j);
            auto it = journal.record(breakpoint, std::views::empty<char>);

            THEN ("the journal contains one element") {
                REQUIRE(journal.size() == 1);
            } AND_THEN ("it points to the recorded deletion") {
                REQUIRE(std::ranges::empty(it->sequence()));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == std::max(i, j));
            }
        }
        WHEN("Recording three deletions at different breakpoints") {

            std::vector breakpoints{{
                libjst::to_breakpoint(source, source.begin() + 3, source.begin() + 4),
                libjst::to_breakpoint(source, source.begin() + 3, source.begin() + 6),
                libjst::to_breakpoint(source, source.begin() + 1, source.begin() + 7)
            }};

            journal.record(breakpoints[0], std::views::empty<char>);
            journal.record(breakpoints[1], std::views::empty<char>);
            journal.record(breakpoints[2], std::views::empty<char>);

            THEN ("the journal contains three element") {
                REQUIRE(journal.size() == 3);
            } AND_THEN ("the elements are ordered by breakpoints in ascending order") {
                auto it = journal.begin();
                REQUIRE(std::ranges::empty(it->sequence()));
                REQUIRE(libjst::low_breakend(*it) == 1);
                REQUIRE(libjst::high_breakend(*it) == 7);
                ++it;
                REQUIRE(std::ranges::empty(it->sequence()));
                REQUIRE(libjst::low_breakend(*it) == 3);
                REQUIRE(libjst::high_breakend(*it) == 4);
                ++it;
                REQUIRE(std::ranges::empty(it->sequence()));
                REQUIRE(libjst::low_breakend(*it) == 3);
                REQUIRE(libjst::high_breakend(*it) == 6);
            }
        }
    }
}

SCENARIO("Recording substitutions", "[basic_multisequence_journal][substitution]") {
    using namespace std::string_literals;

    GIVEN("A basic_multisequence_journal initialized with a std::string source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::basic_multisequence_journal journal{source};

        WHEN("Recording a substitution at different breakpoints") {
            auto i = GENERATE(Catch::Generators::values({0, 8, 16}));
            auto j = GENERATE(Catch::Generators::values({1, 8, 16}));
            std::string alt_sequence('x', std::max(i, j) - i);

            auto breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + j);
            auto it = journal.record(breakpoint, alt_sequence);

            THEN ("the journal contains one element") {
                REQUIRE(journal.size() == 1);
            } AND_THEN ("it points to the recorded deletion") {
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequence));
                REQUIRE(libjst::low_breakend(*it) == i);
                REQUIRE(libjst::high_breakend(*it) == std::max(i, j));
            }
        }
        WHEN("Recording three substitutions at different breakpoints") {

            std::vector breakpoints{{
                libjst::to_breakpoint(source, source.begin() + 3, source.begin() + 4),
                libjst::to_breakpoint(source, source.begin() + 3, source.begin() + 6),
                libjst::to_breakpoint(source, source.begin() + 1, source.begin() + 7)
            }};
            std::vector alt_sequences{ std::string{1, 'x'}, std::string{3, 'x'}, std::string{6, 'x'} };

            journal.record(breakpoints[0], alt_sequences[0]);
            journal.record(breakpoints[1], alt_sequences[1]);
            journal.record(breakpoints[2], alt_sequences[2]);

            THEN ("the journal contains three element") {
                REQUIRE(journal.size() == 3);
            } AND_THEN ("the elements are ordered by breakpoints") {
                auto it = journal.begin();
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequences[2]));
                REQUIRE(libjst::low_breakend(*it) == 1);
                REQUIRE(libjst::high_breakend(*it) == 7);
                ++it;
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequences[0]));
                REQUIRE(libjst::low_breakend(*it) == 3);
                REQUIRE(libjst::high_breakend(*it) == 4);
                ++it;
                REQUIRE(std::ranges::equal(it->sequence(), alt_sequences[1]));
                REQUIRE(libjst::low_breakend(*it) == 3);
                REQUIRE(libjst::high_breakend(*it) == 6);
            }
        }
    }
}

SCENARIO("Fuzzy testing of recording variants", "[basic_multisequence_journal][fuzzy]")
{
    GIVEN("A basic_multisequence_journal initialized with a std::string source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::basic_multisequence_journal journal{source};

        auto i = GENERATE(take(10, Catch::Generators::random(0, 16)));
        auto j = GENERATE(take(10, Catch::Generators::random(0, 16)));
        j = std::max(i, j);
        auto k = GENERATE(take(1, Catch::Generators::random(0, 2))); // 0: substitution, 1: insertion, 2: deletion

        libjst::sequence_breakpoint_t<std::string> breakpoint{};
        std::string alt_sequence{};

        switch (k) {
            case 0: // substitution
                breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + j);
                alt_sequence = std::string(j - i, 'x');
                break;
            case 1: // insertion
                breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + i);
                alt_sequence = std::string(j - i, 'x');
                break;
            case 2: // deletion
                breakpoint = libjst::to_breakpoint(source, source.begin() + i, source.begin() + j);
                alt_sequence = std::string{};
                break;
            default:
                FAIL();
        }

        THEN ("Recording the variant should not throw") {
            REQUIRE_NOTHROW(journal.record(breakpoint, alt_sequence));
        }
    }
}
