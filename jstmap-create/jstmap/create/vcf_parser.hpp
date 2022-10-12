// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>

#include <jstmap/global/jstmap_jst_types.hpp>
// #include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap
{

// std::vector<jst_t> construct_jst_from_vcf(std::filesystem::path const &, std::filesystem::path const &);

// read them all
void construct_jst_from_vcf2(std::filesystem::path const &, std::filesystem::path const &, std::filesystem::path const &);

}  // namespace jstmap
