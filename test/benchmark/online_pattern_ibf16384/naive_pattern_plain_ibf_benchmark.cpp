// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <numeric>

#include <seqan3/test/performance/units.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>

#include <seqan/find.h>

template <typename rcs_store_t>
inline size_t total_bytes(rcs_store_t const & rcs_store) noexcept {
    return std::ranges::size(rcs_store.source()) * rcs_store.size();
}

template <typename ...args_t>
static void naive_pattern(benchmark::State & state, args_t && ...args)
{
    using sequence_t = jstmap::reference_t;
    auto [jst_file, needle_file] = std::tuple{args...};
    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t needle = jstmap::load_queries(needle_file)[0].sequence();

    size_t sequence_count = rcs_store.size();
    size_t batch_size = 16;
    size_t total_runs = sequence_count / batch_size;
    std::vector<sequence_t> batch{batch_size, rcs_store.source()};

    seqan2::Pattern<sequence_t, seqan2::Simple> pattern{needle};
    size_t hit_count{};
    for (auto _ : state)
    {
        hit_count = 0;
        for (size_t run = 0; run < total_runs; ++run) {
            std::ranges::for_each(batch, [&] (auto const & seq) {
                using range_t = std::remove_reference_t<decltype(seq)>;
                seqan2::Finder<range_t> finder{seq};
                while (seqan2::find(finder, pattern)) {
                    ++hit_count;
                }
            });
        }
    }
    state.counters["bytes"] = total_bytes(rcs_store);
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(total_bytes(rcs_store));
    state.counters["#hits"] = hit_count;
}

BENCHMARK_CAPTURE(naive_pattern,
                  online_pattern_plain_needle32,
                  "/home/rahn/workspace/jstmap/build/bench/Release/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle32.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(naive_pattern,
                  online_pattern_plain_needle64,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle64.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(naive_pattern,
                  online_pattern_plain_needle128,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle128.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(naive_pattern,
                  online_pattern_plain_needle256,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle256.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_MAIN();
