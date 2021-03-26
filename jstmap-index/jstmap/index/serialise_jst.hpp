// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map indexer.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/filesystem>

#include <cereal/archives/binary.hpp>

namespace jstmap
{

template <typename jst_t>
void serialise_jst(jst_t const & tree, std::filesystem::path const & output_path)
{
    std::ofstream output_stream{output_path.c_str()};
    cereal::BinaryOutputArchive binary_archive{output_stream};
    tree.save(binary_archive);
}

} // namespace jstmap
