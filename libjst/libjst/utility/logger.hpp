// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides logging utilities for printing more debug information during the JST construction and traversal.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/debug_stream.hpp>

namespace libjst::utility
{

/*!\brief Logs the given arguments to the seqan3::debug_stream.
 * \tparam arguments_t The template parameter pack of the passed arguments; each argument must support the
 *                     formatted output concept to the seqan3::debug_stream.
 *
 * \param[in] line The line number where the log message was issued.
 * \param[in] file_name The name of the file where the log message was issued.
 * \param[in] args The arguments to log to the seqan3::debug_stream.
 *
 * \details
 *
 * Logs a debug message to the seqan3::debug_stream with the file and line number to easily jump to the
 * respective logging message.
 */
template <typename ...arguments_t>
//!\cond
    requires (requires ()
    {
        { seqan3::debug_stream << std::declval<arguments_t>() };
    } && ...)
//!\endcond
constexpr void log_debug(uint32_t const line, const char* file_name, arguments_t && ...args)
{
    seqan3::debug_stream << file_name << ":" << line << ": debug: ";
    (seqan3::debug_stream << ... << std::forward<arguments_t>(args)) << '\n';
}

#if defined(LIBJST_DEBUGGING_ENABLED) && LIBJST_DEBUGGING_ENABLED == 1
    #define LIBJST_LOG_DEBUG(...) libjst::utility::log_debug(__LINE__, __FILE__, __VA_ARGS__)
#else
    #define LIBJST_LOG_DEBUG(...)
#endif // defined(LIBJST_DEBUGGING_ENABLED)

}  // namespace libjst::utility
