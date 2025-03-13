// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <libspm/matcher/shiftor_matcher_restorable.hpp>

#include "fixture_resumable_pattern_ibf.hpp"

namespace just::bench {

BENCHMARK_TEMPLATE_DEFINE_F(fixture_resumable_pattern_ibf, shiftor, capture<&chr22_needle32_ibf2048>)(benchmark::State& state) {
    spm::restorable_shiftor_matcher matcher{needle()};
    run(state, matcher);
}

BENCHMARK_REGISTER_F(fixture_resumable_pattern_ibf, shiftor)
    ->RangeMultiplier(2)->Range(1,1)
    ->UseRealTime();
} // namespace just::bench

BENCHMARK_MAIN();
