// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides sequence tree node label using journaled sequence.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <exception>
#include <type_traits>

#include <libjst/sequence/journaled_sequence.hpp>
#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <std::integral position_t, typename source_t>
        // requires source_t is a viewable range
    class journaled_sequence_label {
    private:

        using journaled_sequence_type = journaled_sequence<source_t>;
        using journaled_sequence_iterator = std::ranges::iterator_t<journaled_sequence_type const>;
        using offset_type = std::make_signed_t<position_t>;

        journaled_sequence_type _journaled_source{}; //!\brief Journal structure managing alternate path sequence.
        position_t _left_position{}; //!\brief Left journaled sequence position marking begin of node label.
        position_t _right_position{}; //!\brief Right journaled sequence position marking end of node label.
        offset_type _offset{}; //!\brief Offset between journaled sequence positions and corresponding source position.

    public:

        using position_type = position_t;
        using size_type = position_t;
        using sequence_type = std::ranges::subrange<journaled_sequence_iterator, journaled_sequence_iterator>;

        static constexpr size_type npos{std::numeric_limits<size_type>::max()};

        journaled_sequence_label() = default;
        journaled_sequence_label(source_t source) noexcept :
            _journaled_source{(source_t &&) source},
            _left_position{0},
            _right_position{static_cast<position_t>(std::ranges::size(_journaled_source))}
        {}

        constexpr sequence_type sequence(size_type const first = 0, size_type const last = npos) const noexcept {
            assert(first <= last);
            size_type max_end = std::min<size_type>(last, std::ranges::size(_journaled_source));
            return std::ranges::subrange{std::ranges::next(_journaled_source.begin(), to_alt_position(first)),
                                         std::ranges::next(_journaled_source.begin(), to_alt_position(max_end))};
        }

        constexpr sequence_type node_sequence() const noexcept {
            return sequence(get_left_position(), get_right_position());
        }

        constexpr sequence_type path_sequence() const noexcept {
            return sequence();
        }

        constexpr position_type get_left_position() const noexcept {
            return _left_position;
        }

        constexpr position_type get_right_position() const noexcept {
            return _right_position;
        }

        constexpr position_type label_size() const noexcept {
            return get_right_position() - get_left_position();
        }

        constexpr void reset_positions(position_type const left_position, position_type const right_position) noexcept {
            assert(left_position <= right_position);
            assert(right_position <= static_cast<position_type>(std::ranges::ssize(path_sequence())));
            _left_position = left_position;
            _right_position = right_position;
        }

        template <typename variant_t>
        constexpr void record_variant(variant_t && variant) {
            record_variant_impl(variant);
            update_label_positions((variant_t &&)variant);
        }

    protected:

        template <typename variant_t>
        constexpr void record_variant_impl(variant_t && variant) {
            auto const alt_position = _journaled_source.begin() + to_alt_position(libjst::low_breakend(variant)); // TODO: replace with low_breakend
            auto alt_seq = libjst::alt_sequence(variant);
            switch (libjst::alt_kind(variant)) {
                case alternate_sequence_kind::replacement: {
                    auto alt_seq_size = std::ranges::size(alt_seq);
                    _journaled_source.replace(alt_position, alt_position + alt_seq_size, std::move(alt_seq));
                    break;
                } case alternate_sequence_kind::deletion: {
                    auto && breakpt = libjst::get_breakpoint(variant);
                    _journaled_source.erase(alt_position, alt_position + libjst::breakpoint_span(breakpt));
                    break;
                } case alternate_sequence_kind::insertion: {
                    _journaled_source.insert(alt_position, std::move(alt_seq));
                    break;
                } default: {
                    //no-op
                }
            }
        }

        constexpr position_type to_alt_position(position_type const ref_position) const noexcept {
            assert(-_offset <= static_cast<offset_type>(ref_position));
            return ref_position + _offset;
        }

        template <typename variant_t>
        constexpr void update_label_positions(variant_t && variant) noexcept {
            position_type const alt_position = to_alt_position(libjst::low_breakend(variant));
            reset_positions(alt_position, alt_position + std::ranges::size(libjst::alt_sequence(variant)));
            _offset += libjst::effective_size(variant);
        }
    };
}  // namespace libjst
