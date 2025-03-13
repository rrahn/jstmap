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

#include <libspm/matcher/shiftor_matcher.hpp>
#include <libjst/rcms/haplotype_viewer.hpp>

template <typename rcs_store_t>
inline size_t total_bytes(rcs_store_t const & rcs_store) noexcept {
    libjst::haplotype_viewer viewer{rcs_store};
    size_t size{};
    for (size_t idx = 0; idx < viewer.size(); ++idx) {
        size += std::ranges::size(viewer[idx]);
    }
    return size;
}

template <typename ...args_t>
static void shiftor_pattern(benchmark::State & state, args_t && ...args)
{
    using sequence_t = jstmap::reference_t;
    auto [jst_file, needle_file] = std::tuple{args...};
    jstmap::rcs_store_t rcs_store = jstmap::load_jst(jst_file);
    sequence_t needle = jstmap::load_queries(needle_file)[0].sequence();

    spm::shiftor_matcher pattern(needle);
    libjst::haplotype_viewer viewer{rcs_store};

    size_t hit_count{};
    for (auto _ : state)
    {
        hit_count = 0;
        for (size_t idx = 0; idx < viewer.size(); ++idx) {
            pattern(viewer[idx], [&] (auto const &) { ++hit_count; });
        }
    }
    state.counters["bytes"] = total_bytes(rcs_store);
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(total_bytes(rcs_store));
}

BENCHMARK_CAPTURE(shiftor_pattern,
                  online_pattern_plain_needle32,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle32.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(shiftor_pattern,
                  online_pattern_plain_needle64,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle64.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(shiftor_pattern,
                  online_pattern_plain_needle128,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle128.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_CAPTURE(shiftor_pattern,
                  online_pattern_plain_needle256,
                  DATADIR"ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.jst",
                  DATADIR"needle256.fa")
                    ->MeasureProcessCPUTime()
                    ->UseRealTime();

BENCHMARK_MAIN();
