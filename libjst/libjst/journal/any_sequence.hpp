// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides polymorphic sequence wrapper.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <exception>
#include <memory>
#include <string_view>

namespace libjst
{

    class bad_sequence_access : public std::exception
    {
    private:

        std::string_view _message{"Bad sequence access. No sequence stored in any_sequence!"};
    public:

        bad_sequence_access() noexcept = default;
        bad_sequence_access(bad_sequence_access const &) noexcept = default;
        bad_sequence_access & operator=(bad_sequence_access const &) noexcept = default;
        virtual ~bad_sequence_access() = default;

        virtual const char* what() const noexcept override
        {
            return _message.data();
        }
    };

    template <typename sequence_t>
    class any_sequence
    {

        /// @name Member types
        /// @{
    public:
        using value_type = sequence_t;

    private:
        class sequence_base
        {
        public:
            constexpr sequence_base() = default;
            virtual ~sequence_base() = default;

            virtual sequence_t get() const = 0;
        };

        template <typename concrete_sequence_t>
            requires std::convertible_to<concrete_sequence_t, sequence_t>
        class sequence_impl final : public sequence_base
        {
        private:
            concrete_sequence_t _sequence{};

        public:
            constexpr sequence_impl() = default;
            constexpr sequence_impl(concrete_sequence_t sequence) : _sequence{std::move(sequence)} {}
            virtual ~sequence_impl() = default;

            constexpr sequence_t get() const override
            {
                return _sequence;
            }
        };
    /// @}

        /// @name Member variables
        /// @{
    private:
        std::unique_ptr<sequence_base> _sequence_ptr{};
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr any_sequence() = default;

        template <typename concrete_sequence_t>
        constexpr any_sequence(concrete_sequence_t sequence)
            noexcept(std::is_nothrow_move_constructible_v<concrete_sequence_t>) :
            _sequence_ptr{std::make_unique<sequence_impl<concrete_sequence_t>>(std::move(sequence))}
        {}
        /// @}

        /// @name Accessor
        /// @{
        constexpr sequence_t value() const
        {
            if (!has_value())
                throw bad_sequence_access{};

            return _sequence_ptr->get();
        }

        constexpr sequence_t operator*() const
        {
            assert(has_value());
            return _sequence_ptr->get();
        }

        constexpr bool has_value() const noexcept
        {
            return _sequence_ptr != nullptr;
        }

        constexpr explicit operator bool() const noexcept
        {
            return has_value();
        }
        /// @}

    };
}  //namespace libjst
