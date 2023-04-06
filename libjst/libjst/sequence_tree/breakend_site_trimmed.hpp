// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides breakend site for sequence graph.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename breakend_site_t>
    class breakend_site_trimmed
    {
    public:

        using delta_reference = typename breakend_site_t::delta_reference;
        using delta_value = typename breakend_site_t::delta_value;
        using index_type = typename breakend_site_t::index_type;
        using value_type = typename breakend_site_t::value_type;
        using position_value_type = libjst::variant_position_t<delta_reference>;

    private:


        breakend_site_t const & _wrappee{};
        position_value_type _max_position{};

    public:

        breakend_site_trimmed() = delete;
        explicit constexpr breakend_site_trimmed(breakend_site_t const & breakend_site,
                                                 position_value_type max_position = std::numeric_limits<position_value_type>::max()) :
            _wrappee{breakend_site},
            _max_position{max_position}
        {
        }

        constexpr delta_reference operator*() const noexcept {
            return *_wrappee;
        }

        constexpr auto get_breakend() const noexcept {
            return _wrappee.get_breakend();
        }

        constexpr bool is_high_end() const noexcept {
            return _wrappee.is_high_end();
        }

        constexpr bool is_low_end() const noexcept {
            return _wrappee.is_low_end();
        }

    private:

        friend constexpr position_value_type
        tag_invoke(std::tag_t<libjst::position>, breakend_site_trimmed const & me) noexcept
        {
            return std::min(libjst::position(me._wrappee), me._max_position);
        }

        constexpr friend bool operator==(breakend_site_trimmed const &, breakend_site_trimmed const &) noexcept = default;
    };

}  // namespace libjst
