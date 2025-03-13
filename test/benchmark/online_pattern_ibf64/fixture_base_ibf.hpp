// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <omp.h>

#include <functional>
#include <ranges>
#include <stack>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/filter_queries.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "fixture_config.hpp"

namespace just::bench
{

    template <auto _config_generator>
    struct capture {
        auto fixture() -> std::invoke_result_t<decltype(*_config_generator)>
        {
            return std::invoke(*_config_generator);
        }
    };

    template <typename capture_t>
    class fixture_base_ibf : public ::benchmark::Fixture,
                             protected capture_t
    {
    protected:
        using sequence_t = jstmap::reference_t;
        using queries_t = std::vector<jstmap::sequence_record_t>;

        size_t processed_bytes;
    private:

        jstmap::rcs_store_t _rcs_store{};
        queries_t _queries{};
        jstmap::search_options _options{};

    public:

        fixture_base_ibf() = default;
        virtual ~fixture_base_ibf() = default;

        virtual void SetUp(::benchmark::State const &) override {
            auto [jst_file, needle_file, ibf_file] = this->fixture();

            _rcs_store = jstmap::load_jst(jst_file);
            _queries = jstmap::load_queries(needle_file);
            _options.index_input_file_path = ibf_file;
            _options.thread_count = 1;
            _options.error_rate = 0.0;
        }

        virtual void TearDown(::benchmark::State & state) override {
            state.counters["bytes"] = processed_bytes;
            state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(state.counters["bytes"]);
        }

        sequence_t const & needle() const noexcept {
            return _queries[0].sequence();
        }

        jstmap::rcs_store_t const & store() const noexcept {
            return _rcs_store;
        }

        auto queries() const noexcept {
            return _queries | std::views::transform([] (auto const & record) -> sequence_t const & {
                return record.sequence();
            });
        }

        constexpr float to_error_rate(int32_t error_count) noexcept {
            _options.error_rate = static_cast<float>(error_count)/100.0 + 0.00001;
            return _options.error_rate;
        }

        inline size_t total_bytes(size_t window_size) noexcept {
            auto tree = store() | libjst::make_volatile()
                                | libjst::labelled()
                                | libjst::coloured()
                                | libjst::trim(window_size - 1)
                                | libjst::prune()
                                | libjst::left_extend(window_size - 1)
                                | libjst::merge();

            return libjst::stats(tree).symbol_count;
        }

        template <typename matcher_t, typename tree_closure_t, typename traverser_factory_t>
        void run(::benchmark::State & state,
                 matcher_t && make_pattern,
                 tree_closure_t && closure,
                 traverser_factory_t && make_traverser)
        {
            _options.thread_count = state.range(0);
            auto queries = make_queries();
            int32_t hit_count{};
            for (auto _ : state)
            {
                auto [bin_size, search_queries] = jstmap::filter_queries(queries, _options);
                auto trees = store() | libjst::chunk(bin_size);

                benchmark::DoNotOptimize(hit_count = execute(trees, make_pattern, closure, search_queries, make_traverser, _options));
                benchmark::ClobberMemory();
            }
        }

        constexpr auto make_queries() const noexcept {
            size_t query_idx{};
            auto queries = _queries
                     | std::views::transform([&] (jstmap::sequence_record_t const & record) {
                        return jstmap::search_query{query_idx++, record};
                     })
                     | seqan3::ranges::to<std::vector>();
            return queries;
        }

        template <typename trees_t, typename matcher_t, typename tree_closure_t, typename search_queries_t, typename traverser_factory_t>
        static int32_t execute(trees_t && trees,
                               matcher_t && make_pattern,
                               tree_closure_t && tree_closure,
                               search_queries_t const & queries,
                               traverser_factory_t && make_traverser,
                               jstmap::search_options const & options) {
            int32_t hit_count = 0;
            std::ptrdiff_t const chunk_count = std::ranges::ssize(trees);

            #pragma omp parallel for num_threads(options.thread_count), shared(trees), firstprivate(make_pattern, make_traverser, tree_closure), schedule(static), reduction(+:hit_count)
            for (std::ptrdiff_t chunk = 0; chunk < chunk_count; ++chunk) {
                if (std::ranges::empty(queries[chunk])) continue;

                auto pattern = make_pattern(queries[chunk] |
                                    std::views::transform([] (jstmap::search_query const & query) {
                                        return std::views::all(query.value().sequence());
                                    }));
                auto tree = trees[chunk] | tree_closure(spm::window_size(pattern));
                auto traverser = make_traverser(tree);
                for (auto it = traverser.begin(); it != traverser.end(); ++it) {
                    auto && cargo = *it;
                    pattern(cargo.sequence(), [&] (auto const &) { ++hit_count; });
                }
            }
            return hit_count;
        }
    };
}  // namespace just::bench

