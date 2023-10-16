// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <libcontrib/matcher/pigeonhole_matcher.hpp>

#include "fixture_oblivious_pattern.hpp"

namespace just::bench {

BENCHMARK_TEMPLATE_DEFINE_F(fixture_oblivious_pattern, pigeonhole, capture<&chr22_sim100x100Ke3>)(benchmark::State& state) {
    run(state, jst::contrib::pigeonhole_matcher(queries(), to_error_rate(state.range(1))));
}

BENCHMARK_REGISTER_F(fixture_oblivious_pattern, pigeonhole)
    ->ArgsProduct({
        benchmark::CreateRange(1, std::thread::hardware_concurrency(), 2),
        {0, 1, 2, 3}
    })
    ->UseRealTime();
} // namespace just::bench

BENCHMARK_MAIN();
