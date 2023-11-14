// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include <libjst/journal/coverage_augmented_breakpoint_multijournal.hpp>

SCENARIO("Initializing a coverage_augmented_breakpoint_multijournal") {
    GIVEN("An empty source") {
        libjst::coverage_augmented_breakpoint_multijournal<std::string> journal{};

        THEN("The journal should be empty") {
            REQUIRE(journal.empty() == true);
            REQUIRE(journal.size() == 0);
        }
    }

    GIVEN("A preinitialised source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::coverage_augmented_breakpoint_multijournal journal{source};

        THEN("The journal should be empty") {
            REQUIRE(journal.empty() == true);
        }

        THEN("The journal' source should be equal to the given source") {
            REQUIRE(journal.source() == source);
        }
    }
}

SCENARIO("Recording a sequence with its coverage") {
    GIVEN("A coverage_augmented_breakpoint_multijournal initialized with a non-empty source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::coverage_augmented_breakpoint_multijournal journal{source};
        using coverage_type = typename decltype(journal)::coverage_type;
        using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

        auto i = GENERATE(0, 4, 16);
        WHEN("I call the record function to insert a sequence and its coverage") {
            std::string sequence{"ACGT"};

            coverage_type coverage{{1, 3}, coverage_domain_type{0, 4}};
            auto breakpoint = libjst::to_breakpoint(journal.source(),
                                                    journal.source().begin() + i,
                                                    journal.source().begin() + i);
            auto it = journal.record(breakpoint, sequence, coverage);

            THEN("The sequence should be recorded in the journal") {
                REQUIRE(journal.empty() == false);
                REQUIRE(journal.size() == 1);
                REQUIRE(std::ranges::equal((*it).sequence(), sequence));
            }

            THEN("The coverage should be associated with the sequence in the journal") {
                REQUIRE((*it).coverage() == coverage);
            }
        }
    }
}

SCENARIO("Recording a deletion with its coverage") {
    GIVEN("A coverage_augmented_breakpoint_multijournal initialized with a non-empty source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::coverage_augmented_breakpoint_multijournal journal{source};
        using coverage_type = typename decltype(journal)::coverage_type;
        using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

        auto i = GENERATE(0, 4, 16);
        WHEN("I call the record function to insert a deletion and its coverage") {
            std::string sequence{""};

            coverage_type coverage{{0, 1}, coverage_domain_type{0, 4}};
            auto breakpoint = libjst::to_breakpoint(journal.source(),
                                                    journal.source().begin() + i,
                                                    journal.source().begin() + std::max(i + 4, 16));
            auto it = journal.record(breakpoint, sequence, coverage);

            THEN("The deletion should be recorded in the journal") {
                REQUIRE(journal.empty() == false);
                REQUIRE(journal.size() == 1);
                REQUIRE(std::ranges::equal((*it).sequence(), sequence));
            }

            THEN("The coverage should be associated with the deletion in the journal") {
                REQUIRE((*it).coverage() == coverage);
            }
        }
    }
}

SCENARIO("Successfully recording various sequence modifications with their coverages") {
    using namespace std::literals;

    GIVEN("A coverage_augmented_breakpoint_multijournal initialized with a non-empty source") {
        std::string source{"AAAACCCCGGGGTTTT"};
        libjst::coverage_augmented_breakpoint_multijournal journal{source};
        using coverage_type = typename decltype(journal)::coverage_type;
        using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

        std::vector alternates{{
            std::tuple{0, 4,  "ACGT"s, coverage_type{{1, 3}, coverage_domain_type{0, 4}}},
            std::tuple{2, 12, ""s, coverage_type{{0}, coverage_domain_type{0, 4}}},
            std::tuple{8, 8, "ACGTACGT"s, coverage_type{{2, 3}, coverage_domain_type{0, 4}}},
            std::tuple{12, 13, "T"s, coverage_type{{0, 3}, coverage_domain_type{0, 4}}},
            std::tuple{16, 16, "ACGT"s, coverage_type{{1, 2}, coverage_domain_type{0, 4}}}
        }};

        std::vector shuffled_alternates{alternates};
        std::shuffle(shuffled_alternates.begin(), shuffled_alternates.end(), std::mt19937{Catch::getSeed()});

        WHEN("I call the record function to insert alternate sequences and their coverages") {
            size_t expected_size = 0;
            for (auto const& [i, j, sequence, coverage] : shuffled_alternates) {
                auto breakpoint = libjst::to_breakpoint(journal.source(),
                                                        journal.source().begin() + i,
                                                        journal.source().begin() + j);
                auto it = journal.record(breakpoint, sequence, coverage);

                THEN("The sequence should be recorded in the journal") {
                    REQUIRE(journal.empty() == false);
                    REQUIRE(journal.size() == ++expected_size);
                    REQUIRE(std::ranges::equal((*it).sequence(), sequence));
                }

                THEN("The coverage should be associated with the sequence in the journal") {
                    REQUIRE((*it).coverage() == coverage);
                }
            }

            THEN("The sequences should be recorded in the journal according to their breakpoints in ascending order") {
                auto it = journal.begin();
                for (auto const& [i, j, sequence, coverage] : alternates) {
                    auto breakpoint = libjst::to_breakpoint(journal.source(),
                                                            journal.source().begin() + i,
                                                            journal.source().begin() + j);
                    REQUIRE(std::ranges::equal((*it).sequence(), sequence));
                    REQUIRE((*it).coverage() == coverage);
                    REQUIRE(libjst::low_breakend(*it) == libjst::low_breakend(breakpoint));
                    REQUIRE(libjst::high_breakend(*it) == libjst::high_breakend(breakpoint));
                    REQUIRE(libjst::breakend_span(*it) == libjst::breakend_span(breakpoint));
                    ++it;
                }
            }
        }
    }
}

