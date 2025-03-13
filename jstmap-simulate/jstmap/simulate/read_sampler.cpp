// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function for loading in a reference sequence for the simulated alignment.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#include <ranges>

#include <libspm/seqan/alphabet.hpp>

#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/simulate/read_sampler.hpp>
#include <jstmap/global/application_logger.hpp>

namespace jstmap
{

    read_sampler::read_sampler(rcs_store_t const & rcs_store) : _rcs_store{rcs_store}
    {}

    libjst::tree_stats
    read_sampler::compute_tree_stats(size_t const size) const {
        log_debug("Compute tree statistics");
        log_debug("Variant count:", _rcs_store.variants().size());
        return libjst::volatile_tree(_rcs_store) | libjst::labelled()
                                                 | libjst::coloured()
                                                 | libjst::trim(size - 1)
                                                 | libjst::prune()
                                                 | libjst::merge()
                                                 | libjst::stats();
    }

    sampled_read_list_type
    read_sampler::sample_reads(libjst::tree_stats const & stats,  size_t const count,  size_t const size) const {
        assert(size > 0);
        assert(count > 0);

        std::mt19937 random_engine{42};
        size_t const window_size = size - 1;
        auto sample_range = compute_sample_range(stats, size);

        log_debug("Window size:", window_size);
        log_debug("Max sample range", sample_range);
        std::vector<size_t> sample_positions = generate_sample_positions(count, sample_range, random_engine);

        sampled_read_list_type sampled_reads{};
        sampled_reads.reserve(sample_positions.size());

        auto sample_tree = libjst::volatile_tree(_rcs_store) | libjst::labelled()
                                                             | libjst::coloured()
                                                             | libjst::trim(window_size)
                                                             | libjst::prune()
                                                             | libjst::left_extend(window_size)
                                                             | libjst::merge()
                                                             | libjst::seek();

        libjst::tree_traverser_base path{sample_tree};
        auto next_sample_position = sample_positions.begin();
        auto tree_it = path.begin();
        size_t last_position = window_size;
        while (tree_it != path.end()) {
            auto cargo = *tree_it;
            auto label = cargo.sequence();
            assert(std::ranges::size(label) >= window_size);
            size_t current_position = last_position + std::ranges::size(label) - window_size;
            if (*next_sample_position <= current_position) { // found sample position
                log_debug("====== Found sample position ======");
                log_debug("Current position:", current_position);
                log_debug("Sample position:", *next_sample_position);
                log_debug("Label size extended:", std::ranges::size(label));
                log_debug("Label size normal:", std::ranges::size(label) - window_size);
                log_debug("Sample id:", std::ranges::distance(sample_positions.begin(), next_sample_position));

                std::ptrdiff_t sample_end_offset = current_position - *next_sample_position;
                std::ptrdiff_t sample_end = std::ranges::size(label) - sample_end_offset;
                std::ptrdiff_t sample_begin = sample_end - size;

                log_debug("Sample end offset:", sample_end_offset);
                log_debug("Sample end:", sample_end);
                log_debug("Sample begin:", sample_begin);

                assert(sample_end_offset < (std::ranges::size(label) - window_size));
                assert(sample_end >= size);

                auto read_begin = std::ranges::next(std::ranges::begin(label), sample_begin);
                auto read_end = std::ranges::next(std::ranges::begin(label), sample_end);
                read_type sample{read_begin, read_end};
                using spm::operator""_dna5;
                if (std::ranges::none_of(sample, [] (auto && symbol) { return symbol == 'N'_dna5; })) {
                    sampled_reads.emplace_back(std::move(sample),
                                               match_position{.tree_position = cargo.position(),
                                                              .label_offset = sample_end_offset});
                    validate_sample(sampled_reads.back());
                }
                // assert();

                if (++next_sample_position == sample_positions.end())
                    break;
            } else {
                ++tree_it;
                last_position = current_position;
            }
        }
        return sampled_reads;
    }

    std::vector<size_t> read_sampler::generate_sample_positions(size_t const count,
                                                                std::pair<size_t, size_t> const max_sample_range,
                                                                std::mt19937 & random_engine) const {

        auto [from, to] = max_sample_range;
        std::uniform_int_distribution<size_t> sample_dist{from, to};
        std::vector<size_t> tmp{};
        tmp.resize(count);
        std::ranges::generate(tmp, [&] () { return sample_dist(random_engine); });
        std::ranges::sort(tmp);
        return tmp;
    }

    std::pair<size_t, size_t> read_sampler::compute_sample_range(libjst::tree_stats const & stats,
                                                                 size_t const sample_size) const noexcept {
        auto const & reference = _rcs_store.source();
        auto first_it = std::ranges::find_if_not(reference, [] (auto const nucleotide) {
            return seqan3::to_char(nucleotide) == 'N';
        });

        auto last_it = std::ranges::find_if_not(reference | std::views::reverse, [] (auto const nucleotide) {
            return seqan3::to_char(nucleotide) == 'N';
        }).base();

        size_t left_overhead = std::ranges::distance(std::ranges::begin(reference), first_it);
        size_t right_overhead = std::ranges::distance(last_it, std::ranges::end(reference));
        log_debug("Max subtree depth:", stats.max_subtree_depth);
        log_debug("Left overhead:", left_overhead);
        log_debug("Right overhead:", right_overhead);
        log_debug("Symbol count:", stats.symbol_count);
        assert(stats.symbol_count > (left_overhead + right_overhead + sample_size));

        return std::pair{left_overhead + sample_size, stats.symbol_count - right_overhead};
    }

    bool read_sampler::validate_sample(sampled_read_type const & sample) const noexcept {
        auto && [read, match_position] = sample;
        size_t sample_size = std::ranges::size(read);
        size_t window_size = sample_size - 1;
        auto validation_tree = libjst::volatile_tree(_rcs_store) | libjst::labelled()
                                                                 | libjst::coloured()
                                                                 | libjst::trim(window_size)
                                                                 | libjst::prune()
                                                                 | libjst::left_extend(window_size)
                                                                 | libjst::merge()
                                                                 | libjst::seek();

        auto node = validation_tree.seek(match_position.tree_position);
        auto cargo = *node;
        auto label = cargo.sequence();
        size_t end_position = std::ranges::size(label) - match_position.label_offset;
        size_t begin_position = end_position - sample_size;
        auto label_begin = std::ranges::next(std::ranges::begin(label), begin_position);
        auto label_end = std::ranges::next(std::ranges::begin(label), end_position);
        bool test = std::ranges::equal(read, std::ranges::subrange{label_begin, label_end});
        if (!test)
            log_warn("The sampled read does not corresponds with the sequence at the given seek position");
        return test;
    }
} // namespace jstmap
