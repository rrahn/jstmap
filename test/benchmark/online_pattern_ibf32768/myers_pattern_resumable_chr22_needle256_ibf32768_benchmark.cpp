// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <libspm/matcher/myers_matcher_restorable.hpp>

#include "fixture_resumable_pattern_ibf.hpp"

namespace just::bench {

BENCHMARK_TEMPLATE_DEFINE_F(fixture_resumable_pattern_ibf, myers, capture<&chr22_needle256_ibf32768>)(benchmark::State& state) {
    run(state, spm::restorable_myers_matcher(needle(), (size_t) state.range(1)));
}

BENCHMARK_REGISTER_F(fixture_resumable_pattern_ibf, myers)
    ->ArgsProduct({
        benchmark::CreateRange(1, 1, /*multi=*/2),
        {0, 1, 2, 3}
    })
    ->UseRealTime();
} // namespace just::bench

BENCHMARK_MAIN();
