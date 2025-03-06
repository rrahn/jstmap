// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for sequence variants.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <ranges>
#include <type_traits>

#include <libjst/utility/tag_invoke.hpp>

#include <libjst/variant/alternate_sequence_kind.hpp>
#include <libjst/reference_sequence/sequence_breakpoint_concept.hpp>

namespace libjst
{
    // TODO remove me!
    struct breakpoint;

    enum struct breakpoint_end : uint8_t {
        right = 0,
        left = 1,
        low,
        high
    };

    // ----------------------------------------------------------------------------
    // Operation CPOs for variants
    // ----------------------------------------------------------------------------

    // position
    namespace _variant_position {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t && var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } position;
    } // namespace _variant_position
    using _variant_position::position;

    template <typename variant_t>
    using variant_position_t = std::invoke_result_t<_variant_position::_cpo, variant_t>;

    // insertion
    namespace _variant_insertion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } insertion;
    } // namespace insertion
    using _variant_insertion::insertion;

    template <typename variant_t>
    using variant_insertion_t = std::invoke_result_t<_variant_insertion::_cpo, variant_t>;

    // deletion
    namespace _variant_deletion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } deletion;
    }
    using _variant_deletion::deletion;

    template <typename variant_t>
    using variant_deletion_t = std::invoke_result_t<_variant_deletion::_cpo, variant_t>;

    //is_deletion
    namespace _is_deletion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return libjst::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!libjst::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (libjst::tag_invocable<libjst::tag_t<libjst::insertion>, variant_t const &> &&
                              libjst::tag_invocable<libjst::tag_t<libjst::deletion>, variant_t const &>)
                    return libjst::deletion(var) > 0 && std::ranges::size(libjst::insertion(var)) == 0;
                else
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_deletion(variant) found!");
            }

        } is_deletion;
    }
    using _is_deletion::is_deletion;

    //is_insertion
    namespace _is_insertion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return libjst::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!libjst::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (libjst::tag_invocable<libjst::tag_t<libjst::insertion>, variant_t const &> &&
                              libjst::tag_invocable<libjst::tag_t<libjst::deletion>, variant_t const &>)
                    return libjst::deletion(var) == 0 && std::ranges::size(libjst::insertion(var)) > 0;
                else
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_insertion(variant) found!");
            }

        } is_insertion;
    }
    using _is_insertion::is_insertion;

    //is_replacement
    namespace _is_replacement {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return libjst::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!libjst::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (libjst::tag_invocable<libjst::tag_t<libjst::insertion>, variant_t const &> &&
                              libjst::tag_invocable<libjst::tag_t<libjst::deletion>, variant_t const &>)
                {
                    return libjst::deletion(var) == std::ranges::size(libjst::insertion(var)) &&
                           libjst::deletion(var) > 0;
                } else {
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_replacement(variant) found!");
                }
            }

        } is_replacement;
    }
    using _is_replacement::is_replacement;

    // coverage
    namespace _variant_coverage {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } coverage;
    }
    using _variant_coverage::coverage;

    template <typename variant_t>
    using variant_coverage_t = std::remove_cvref_t<std::invoke_result_t<_variant_coverage::_cpo, variant_t>>;

    // ----------------------------------------------------------------------------
    // CPO libjst::get_breakpoint
    // ----------------------------------------------------------------------------

    // New version with updated name!
    namespace _get_breakpoint {
        inline constexpr struct _cpo  {
            template <typename breakpoint_t>
                requires libjst::tag_invocable<_cpo, breakpoint_t>
            constexpr auto operator()(breakpoint_t &&bp) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, breakpoint_t>)
                -> libjst::tag_invoke_result_t<_cpo, breakpoint_t>
            {
                return libjst::tag_invoke(_cpo{}, (breakpoint_t &&)bp);
            }
        } get_breakpoint;
    } // namespace _get_breakpoint
    using _get_breakpoint::get_breakpoint;


    template <typename breakpoint_t>
    using breakend_t =
        std::remove_cvref_t<
            std::common_type_t<low_breakend_t<breakpoint_t>, high_breakend_t<breakpoint_t>>
        >;

    // ----------------------------------------------------------------------------
    // CPO libjst::breakpoint_span
    // ----------------------------------------------------------------------------
    namespace _breakpoint_span {
        inline constexpr struct _cpo  {
            template <typename sequence_alternative_t>
                requires libjst::tag_invocable<_cpo, sequence_alternative_t>
            constexpr auto operator()(sequence_alternative_t &&alt) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, sequence_alternative_t>)
                -> libjst::tag_invoke_result_t<_cpo, sequence_alternative_t>
            {
                return libjst::tag_invoke(_cpo{}, (sequence_alternative_t &&)alt);
            }

        private:

            template <typename breakpoint_t>
                requires libjst::tag_invocable<libjst::tag_t<libjst::low_breakend>, breakpoint_t> &&
                         libjst::tag_invocable<libjst::tag_t<libjst::high_breakend>, breakpoint_t>
            constexpr friend auto tag_invoke(_cpo, breakpoint_t const & bp)
                noexcept(libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::low_breakend>, breakpoint_t const &> &&
                         libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::high_breakend>, breakpoint_t const &>)
                -> std::common_type_t<libjst::tag_invoke_result_t<libjst::tag_t<libjst::high_breakend>, breakpoint_t const &>,
                                      libjst::tag_invoke_result_t<libjst::tag_t<libjst::low_breakend>, breakpoint_t const &>>
            {
                assert(libjst::low_breakend(bp) <= libjst::high_breakend(bp));
                return libjst::high_breakend(bp) - libjst::low_breakend(bp);
            }
        } breakpoint_span;
    } // namespace _breakpoint_span
    using _breakpoint_span::breakpoint_span;

    template <typename sequence_alternative_t>
    using breakpoint_span_t = std::remove_cvref_t<std::invoke_result_t<_breakpoint_span::_cpo, sequence_alternative_t>>;

    // ----------------------------------------------------------------------------
    // CPO libjst::left_breakpoint
    // ----------------------------------------------------------------------------
    namespace _left_breakpoint {
        inline constexpr struct _cpo  {
            template <typename sequence_variant_t>
                requires libjst::tag_invocable<_cpo, sequence_variant_t>
            constexpr auto operator()(sequence_variant_t &&variant) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, sequence_variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, sequence_variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (sequence_variant_t &&)variant);
            }
        } left_breakpoint;
    } // namespace _left_breakpoint
    using _left_breakpoint::left_breakpoint;

    template <typename sequence_variant_t>
    using left_breakpoint_t =
        std::remove_cvref_t<std::invoke_result_t<_left_breakpoint::_cpo, sequence_variant_t>>;

    // ----------------------------------------------------------------------------
    // CPO libjst::right_breakpoint
    // ----------------------------------------------------------------------------
    namespace _right_breakpoint {
        inline constexpr struct _cpo  {
            template <typename sequence_variant_t>
                requires libjst::tag_invocable<_cpo, sequence_variant_t>
            constexpr auto operator()(sequence_variant_t &&variant) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, sequence_variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, sequence_variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (sequence_variant_t &&)variant);
            }
        private:

            template <typename sequence_variant_t>
                requires libjst::tag_invocable<libjst::tag_t<libjst::left_breakpoint>, sequence_variant_t> &&
                         libjst::tag_invocable<libjst::tag_t<libjst::breakpoint_span>, sequence_variant_t>
            constexpr friend auto tag_invoke(_cpo, sequence_variant_t &&variant)
                noexcept(libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::left_breakpoint>, sequence_variant_t> &&
                         libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::breakpoint_span>, sequence_variant_t>)
                -> libjst::tag_invoke_result_t<libjst::tag_t<libjst::left_breakpoint>, sequence_variant_t>
            {
                using breakpoint_t = libjst::tag_invoke_result_t<libjst::tag_t<libjst::left_breakpoint>, sequence_variant_t>;
                using value_t = typename breakpoint_t::value_type;
                return breakpoint_t{libjst::left_breakpoint(variant).value() +
                                    static_cast<value_t>(libjst::breakpoint_span(variant)),
                                    breakpoint_end::right};
            }
        } right_breakpoint;
    } // namespace _right_breakpoint
    using _right_breakpoint::right_breakpoint;

    template <typename sequence_variant_t>
    using variant_breakpoint_t = std::remove_cvref_t<
                std::common_type_t<std::invoke_result_t<_left_breakpoint::_cpo, sequence_variant_t>,
                                   std::invoke_result_t<_right_breakpoint::_cpo, sequence_variant_t>>>;


    // ----------------------------------------------------------------------------
    // CPO libjst::alt_sequence
    // ----------------------------------------------------------------------------
    namespace _alt_sequence {
        inline constexpr struct _cpo  {
            template <typename sequence_alternative_t>
                requires libjst::tag_invocable<_cpo, sequence_alternative_t>
            constexpr auto operator()(sequence_alternative_t &&alt) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, sequence_alternative_t>)
                -> libjst::tag_invoke_result_t<_cpo, sequence_alternative_t>
            {
                return libjst::tag_invoke(_cpo{}, (sequence_alternative_t &&)alt);
            }
        } alt_sequence;
    } // namespace _alt_sequence
    using _alt_sequence::alt_sequence;

    template <typename sequence_alternative_t>
    using alt_sequence_t = std::invoke_result_t<_alt_sequence::_cpo, sequence_alternative_t>;

    // ----------------------------------------------------------------------------
    // CPO libjst::effective_size
    // ----------------------------------------------------------------------------
    namespace _effective_size {
        inline constexpr struct _cpo  {
            // If tag_invocable
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }

        private:
            template <typename variant_t>
                requires libjst::tag_invocable<libjst::tag_t<libjst::breakpoint_span>, variant_t> &&
                         libjst::tag_invocable<libjst::tag_t<libjst::alt_sequence>, variant_t>
            constexpr friend auto tag_invoke(_cpo, variant_t && var)
                noexcept(libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::breakpoint_span>, variant_t> &&
                         libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::alt_sequence>, variant_t>)
                -> std::common_type_t<libjst::tag_invoke_result_t<libjst::tag_t<libjst::breakpoint_span>, variant_t>,
                                      std::ranges::range_difference_t<
                                            libjst::tag_invoke_result_t<libjst::tag_t<libjst::alt_sequence>, variant_t>>>
            {
                return std::ranges::distance(libjst::alt_sequence(var)) - libjst::breakpoint_span(var);
            }
        } effective_size;
    } // namespace _effective_size
    using _effective_size::effective_size;

    // ----------------------------------------------------------------------------
    // CPO libjst::alt_kind
    // ----------------------------------------------------------------------------
    namespace _alt_kind {
        inline constexpr struct _cpo  {
            // If tag_invocable
            template <typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_t>
            constexpr alternate_sequence_kind operator()(variant_t &&var) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_t>)
            {
                return libjst::tag_invoke(_cpo{}, (variant_t &&)var);
            }

        private:
            template <typename variant_t>
                requires libjst::tag_invocable<libjst::tag_t<libjst::effective_size>, variant_t>
            constexpr friend auto tag_invoke(_cpo, variant_t && var)
                noexcept(libjst::is_nothrow_tag_invocable_v<libjst::tag_t<libjst::effective_size>, variant_t>)
                -> alternate_sequence_kind
            {
                if (libjst::effective_size(var) < 0)
                    return alternate_sequence_kind::deletion;
                else if (libjst::effective_size(var) == 0)
                    return alternate_sequence_kind::replacement;
                else
                    return alternate_sequence_kind::insertion;
            }
        } alt_kind;
    } // namespace _alt_kind
    using _alt_kind::alt_kind;

    // ----------------------------------------------------------------------------
    // Operation CPOs for variant stores
    // ----------------------------------------------------------------------------

    namespace _insert {
        inline constexpr struct _cpo  {
            template <typename variant_store_t, typename variant_t>
                requires libjst::tag_invocable<_cpo, variant_store_t &, variant_t>
            constexpr auto operator()(variant_store_t & store, variant_t && variant) const
                noexcept(libjst::is_nothrow_tag_invocable_v<_cpo, variant_store_t &, variant_t>)
                -> libjst::tag_invoke_result_t<_cpo, variant_store_t &, variant_t>
            {
                return libjst::tag_invoke(_cpo{}, store, (variant_t &&)variant);
            }
        } insert;
    }
    using _insert::insert;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename variant_t>
    concept sequence_variant = requires
    (variant_t const & variant)
    {
        { libjst::position(variant) } -> std::integral;
        { libjst::deletion(variant) } -> std::integral;
        { libjst::insertion(variant) } -> std::ranges::forward_range;
        { libjst::is_replacement(variant) } -> std::same_as<bool>;
        { libjst::is_insertion(variant) } -> std::same_as<bool>;
        { libjst::is_deletion(variant) } -> std::same_as<bool>;
    };

    template <typename variant_t>
    concept covered_sequence_variant = sequence_variant<variant_t> && requires
    (variant_t const & variant)
    {
        { libjst::coverage(variant) } -> std::ranges::random_access_range;
    };

    template <typename store_t>
    concept sequence_variant_store = std::ranges::random_access_range<store_t>
    && requires
    (store_t const & store, std::ranges::range_value_t<store_t> const & variant)
    {
        requires sequence_variant<std::ranges::range_value_t<store_t>>;
        requires sequence_variant<std::ranges::range_reference_t<store_t>>;

        // { store.insert(variant) } -> std::random_access_iterator;
    }
    ;

    template <typename store_t>
    concept covered_sequence_variant_store = sequence_variant_store<store_t>
    && requires
    {
        requires covered_sequence_variant<std::ranges::range_value_t<store_t>>;
        requires covered_sequence_variant<std::ranges::range_reference_t<store_t>>;
    }
    ;

}  // namespace libjst
