// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string_view>
#include <vector>

#include <seqan3/alphabet/gap/gapped.hpp>

namespace libjst::test
{

static std::vector<seqan3::gapped<char>> make_gapped(std::string_view const seq)
{
    std::vector<seqan3::gapped<char>> tmp{};
    tmp.reserve(seq.size());

    std::for_each(seq.begin(), seq.end(), [&] (char const c)
    {
        if (c == '-')
            tmp.emplace_back(seqan3::gap{});
        else
            tmp.emplace_back(c);
    });

    return tmp;
}
} // namespace libjst
