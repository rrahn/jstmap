// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <filesystem>
#include <vector>

#include <jstmap/global/jstmap_types.hpp>

namespace jstmap
{

inline sequence_collection_t load_sequences(std::filesystem::path const &);

} // namespace jstmap
