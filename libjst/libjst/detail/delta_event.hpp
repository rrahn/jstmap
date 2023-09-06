// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_event.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <cereal/types/variant.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/detail/delta_kind_deletion.hpp>
#include <libjst/detail/delta_kind_insertion.hpp>
#include <libjst/detail/delta_kind_snp.hpp>
#include <libjst/detail/delta_kind_substitution.hpp>
#include <libjst/reference_position.hpp>

namespace libjst::detail
{

/*!\brief A delta event represents a difference between a target sequence and a reference sequence.
 *
 * \tparam alphabet_t The alphabet_type used to store inserted sequence information; must model
 *                    seqan3::seqan3::semialphabet.
 *
 * \details
 *
 * A delta event represents a single difference between a target sequence and a reference sequence.
 * In a referentially compressed sequence the target sequence is decomposed into a collection of such delta events,
 * that represent all differences between itself and the respective reference sequence.
 */
template <seqan3::semialphabet alphabet_t>
class delta_event
{
public:
    /*!\name Associated types
     * \{
     */
    using substitution_type = delta_kind_substitution<alphabet_t>; //!< The type of the substitution.
    using snp_type = delta_kind_snp<alphabet_t>; //!< The type of the substitution.
    using insertion_type = delta_kind_insertion<alphabet_t>; //!< The type of the insertion.
    using deletion_type = delta_kind_deletion; //!< The type of the deletion.

    using alphabet_type = alphabet_t; //!< The alphabet type.
    using segment_type = std::span<alphabet_type const>; //!< The segment type.
    using size_type = size_t; //!< The size type.
    using position_type = reference_position; //!< The position type.
    //!\brief The variant type over the three different delta kinds.
    using delta_variant_type = std::variant<insertion_type, snp_type, substitution_type, deletion_type>;
    //!\}

private:
    position_type _position{}; //!< The position of the event.
    delta_variant_type _delta_variant{}; //!< The variant holding one of the valid event types.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    delta_event() = default; //!< Default.
    delta_event(delta_event const &) = default; //!< Default.
    delta_event(delta_event &&) = default; //!< Default.
    delta_event & operator=(delta_event const &) = default; //!< Default.
    delta_event & operator=(delta_event &&) = default; //!< Default.
    ~delta_event() = default; //!< Default.

    /*!\brief Constructs a new delta event from a position and a delta event kind.
     *
     * \param[in] position The position of the delta event.
     * \param[in] kind The kind of the delta event.
     */
    explicit delta_event(position_type position, delta_variant_type kind) :
        _position{std::move(position)},
        _delta_variant{std::move(kind)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the delta event position.
    constexpr position_type position() const noexcept
    {
        return _position;
    }

    //!\brief Returns the delta event variant.
    constexpr delta_variant_type const & delta_variant() const noexcept
    {
        return _delta_variant;
    }
    //!\}

    /*!\name Event operations
     * \{
     */
    //!\brief Returns `true` if this event is a libjst::detail::delta_kind_deletion, `false` otherwise.
    bool is_deletion() const noexcept
    {
        return holds_delta_kind<deletion_type>();
    }

    //!\brief Returns `true` if this event is an libjst::detail::delta_kind_insertion, `false` otherwise.
    bool is_insertion() const noexcept
    {
        return holds_delta_kind<insertion_type>();
    }

    //!\brief Returns `true` if this event is an libjst::detail::delta_kind_substitution, `false` otherwise.
    bool is_substitution() const noexcept
    {
        return is_snp() || holds_delta_kind<substitution_type>();
    }

    bool is_snp() const noexcept
    {
        return holds_delta_kind<snp_type>();
    }

    /*!\brief Returns the deletion size of this event.
     *
     * \returns The size of the deletion.
     *
     * \details
     *
     * The deletion size corresponds to either the length of the deletion or the length of
     * the substituted sequence. If this event is an insertion, the deletion size is 0.
     */
    constexpr size_type deletion_size() const noexcept
    {
        return std::visit([] (auto & event_kind) -> size_type
        {
            return seqan3::detail::multi_invocable
            {
                /*case:*/   [] (substitution_type const & e) { return e.value().size(); },
                /*case:*/   [] (deletion_type const & e) { return e.value(); },
                /*case:*/   [] (snp_type const &) { return 1; },
                /*default:*/[] (...) { return 0; }
            }(event_kind);
        }, delta_variant());
    }

    /*!\brief Returns the insertion size of this event.
     *
     * \returns The size of the insertion.
     *
     * \details
     *
     * The insertion size corresponds to either the length of the insertion or the length of
     * the substituted sequence. If this event is a deletion, the insertion size is 0.
     */
    constexpr size_type insertion_size() const noexcept
    {
        return std::visit([] (auto & event_kind) -> size_type
        {
            return seqan3::detail::multi_invocable
            {
                /*case:*/   [] (substitution_type const & e) { return e.value().size(); },
                /*case:*/   [] (insertion_type const & e) { return e.value().size(); },
                /*case:*/   [] (snp_type const &) { return 1; },
                /*default:*/[] (...) { return 0; }
            }(event_kind);
        }, delta_variant());
    }

    /*!\brief Returns the associated event sequence.
     *
     * \returns The sequence associated with this event.
     *
     * \details
     *
     * The associated sequence is only returned for the insertion and the substitution and is empty for the deletion.
     */
    constexpr segment_type sequence() const noexcept
    {
        return std::visit([] (auto const & event_kind) -> segment_type
        {
            return seqan3::detail::multi_invocable
            {
                /*case:*/   [] (deletion_type const &) { return segment_type{}; },
                /*default:*/[] (auto const & e) { return segment_type{e.value()}; }
            }(event_kind);
        }, delta_variant());
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    //!\brief Compares two delta events for equality.
    bool operator==(delta_event const &) const = default;

    //!\brief Compares the positions of two delta events and returns a std::weak_ordering.
    std::weak_ordering operator<=>(delta_event const & rhs) const noexcept
    {
        return position() <=> rhs.position();
    }
    //!\}

    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves this delta event to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in, out] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_position.offset, _delta_variant);
    }

    /*!\brief Loads this delta event from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in, out] archive The archive to serialise this object from.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(_position.offset, _delta_variant);
    }
    //!\}

private:
    /*!\brief Tests wether the delta variant holds the given delta kind at the moment.
     *
     * \tparam delta_kind_t The type of the delta kind.
     *
     * \returns `true` if the delta variant holds the queried delta kind, `false` otherwise.
     */
    template <typename delta_kind_t>
    constexpr bool holds_delta_kind() const noexcept
    {
        return std::holds_alternative<delta_kind_t>(delta_variant());
    }
};

/*!\brief Formatted output operator for the libjst::detail::delta_event.
 * \relates libjst::detail::delta_event
 *
 * \param[in] stream The basic ostream.
 * \param[in] event The delta event to output to the stream.
 *
 * \returns A reference to the given stream.
 */
template <typename char_t, typename char_traits_t, typename alphabet_t>
inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                              delta_event<alphabet_t> const & event)
{
    using namespace std::literals;

    stream << "(" << event.position() << ", "
           << std::visit([&] (auto event) -> std::string
           {
                auto to_string = [] <std::ranges::input_range range_t> (range_t && range) -> std::string
                {
                    return range | seqan3::views::to_char | seqan3::views::to<std::string>;
                };

                return seqan3::detail::multi_invocable
                {
                    [&] (delta_kind_substitution<alphabet_t> e) -> std::string
                    {
                        return  "sub: " + to_string(e.value());
                    },
                    [&] (delta_kind_snp<alphabet_t> e) -> std::string
                    {
                        return  "snp: " + to_string(e.value());
                    },
                    [&] (delta_kind_insertion<alphabet_t> e) -> std::string
                    {
                        return "ins: " + to_string(e.value());
                    },
                    [&] (delta_kind_deletion e) { return "del: " + std::to_string(e.value()); },
                    [&] (...) { return "N/A"; }
                }(event);
           }, event.delta_variant())
           << ")";
    return stream;
}

} // namespace libjst::detail
