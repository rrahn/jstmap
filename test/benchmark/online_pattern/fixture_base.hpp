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

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/stats.hpp>
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
    class fixture_base : public ::benchmark::Fixture,
                         protected capture_t
    {
    protected:
        using sequence_t = jstmap::reference_t;

        size_t processed_bytes;
    private:

        jstmap::rcs_store_t _rcs_store{};
        sequence_t _needle{};

    public:

        fixture_base() = default;
        virtual ~fixture_base() = default;

        virtual void SetUp(::benchmark::State const &) override {
            auto [jst_file, needle_file] = this->fixture();

            _rcs_store = jstmap::load_jst(jst_file);
            _needle = jstmap::load_queries(needle_file)[0].sequence();
        }

        virtual void TearDown(::benchmark::State & state) override {
            state.counters["bytes"] = processed_bytes;
            state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(state.counters["bytes"]);
        }

        sequence_t const & needle() const noexcept {
            return _needle;
        }

        jstmap::rcs_store_t const & store() const noexcept {
            return _rcs_store;
        }

        template <typename tree_closure_t>
        inline size_t total_bytes(tree_closure_t && tree_closure) noexcept {
            return libjst::stats(tree_closure(store() | libjst::make_volatile())).symbol_count;
        }

        template <typename matcher_t, typename tree_closure_t, typename traverser_factory_t>
        void run(::benchmark::State & state,
                 matcher_t & matcher,
                 tree_closure_t && closure,
                 traverser_factory_t && make_traverser)
        {
            auto trees = libjst::chunk(store(), get_chunk_size(state.range(0)))
                       | std::views::transform([&] (auto && tree) { return closure(tree); });

            int32_t hit_count{};
            for (auto _ : state)
            {
                benchmark::DoNotOptimize(hit_count = execute(trees, matcher, (traverser_factory_t &&) make_traverser));
                benchmark::ClobberMemory();
            }
        }

        constexpr size_t get_chunk_size(size_t const thread_count) noexcept {
            return (std::ranges::size(store().source()) + thread_count - 1) / thread_count;
        }

        template <typename trees_t, typename matcher_t, typename traverser_factory_t>
        static int32_t execute(trees_t && trees, matcher_t & matcher, traverser_factory_t && make_traverser) {
            int32_t hit_count = 0;
            std::ptrdiff_t const chunk_count = std::ranges::ssize(trees);

            #pragma omp parallel for shared(trees), firstprivate(matcher, make_traverser), schedule(static), reduction(+:hit_count)
            for (std::ptrdiff_t chunk = 0; chunk < chunk_count; ++chunk) {
                auto tree = trees[chunk];
                auto traverser = make_traverser(tree);
                for (auto it = traverser.begin(); it != traverser.end(); ++it) {
                    auto && cargo = *it;
                    matcher(cargo.sequence(), [&] (auto const &) { ++hit_count; });
                }
            }
            return hit_count;
        }
    };
}  // namespace just::bench

