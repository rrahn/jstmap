// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map creator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace seqan3
{

class argument_parser;

} // namespace seqan3

namespace jstmap
{

int create_main(seqan3::argument_parser &);

} // namespace jstmap
