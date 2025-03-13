// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides globally defined type aliases for the jstmap tools.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

#include <libspm/seqan/alphabet.hpp>

namespace jstmap
{

//!\brief The sequence type loaded from the disk.
using raw_sequence_t = std::vector<spm::dna5>;

template <typename char_t, seqan3::arithmetic score_t = int8_t>
class scoring_scheme : public seqan3::scoring_scheme_base<scoring_scheme<char_t, score_t>, char_t, score_t>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = seqan3::scoring_scheme_base<scoring_scheme, char_t, score_t>;

    //!\brief Befriend base_t so it can access itself through this derived type.
    friend base_t;

public:
    //!\privatesection
    //!\copydoc scoring_scheme_base::matrix_type
    using typename base_t::matrix_type;
    //!\publicsection

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\copydoc scoring_scheme_base::scoring_scheme_base()
    constexpr scoring_scheme() noexcept = default;
    //!\copydoc scoring_scheme_base::scoring_scheme_base(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    template <seqan3::arithmetic score_arg_t>
    constexpr scoring_scheme(seqan3::match_score<score_arg_t> const ms, seqan3::mismatch_score<score_arg_t> const mms)
        : base_t{ms, mms}
    {}
    //!\copydoc scoring_scheme_base::scoring_scheme_base(matrix_type const & matrix)
    constexpr scoring_scheme(matrix_type const & matrix) noexcept
        : base_t{matrix}
    {}
    //!\}
};

//!\brief Default constructed objects deduce to `int8_t`.
scoring_scheme() -> scoring_scheme<spm::dna5, int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`.
 * To use a larger type, specify the template argument manually.
 */
template <seqan3::arithmetic score_arg_type>
scoring_scheme(seqan3::match_score<score_arg_type>,
               seqan3::mismatch_score<score_arg_type>) -> scoring_scheme<spm::dna5, int8_t>;


struct sequence_input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = spm::dna5;
    using sequence_legal_alphabet = spm::dna15; // conversion?
};
}  // namespace jstmap
