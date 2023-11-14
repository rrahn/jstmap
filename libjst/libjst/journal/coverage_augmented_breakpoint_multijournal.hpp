// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides basic pansequence journal.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/journal/breakpoint_multijournal.hpp>
#include <libjst/coverage/int_coverage.hpp>

namespace libjst
{

    // TODO: add concept for coverage
    template <libjst::reference_sequence source_t>
    class coverage_augmented_breakpoint_multijournal
    {
        /// @name Member types
        /// @{
    private:
        template <typename, typename>
        class record_impl;

        template <bool>
        class iterator_impl;

        using base_journal_t = breakpoint_multijournal<source_t>;
    public:
        using source_type = typename base_journal_t::source_type;
        using sequence_type = typename base_journal_t::sequence_type;
        using breakpoint_type = typename base_journal_t::breakpoint_type;
        using coverage_type = int_coverage<uint32_t>;
        using iterator = iterator_impl<true>;

    private:
        using coverages_type = std::vector<coverage_type>;
        ///@}

        /// @name Member variables
        /// @{
    private:

        base_journal_t _journal{};
        coverages_type _coverages{};
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr coverage_augmented_breakpoint_multijournal() = default;

        constexpr explicit coverage_augmented_breakpoint_multijournal(source_type source)
            noexcept(std::is_nothrow_constructible_v<base_journal_t, source_type>)
            : _journal{std::move(source)}
        {
        }

        constexpr source_type const & source() const & noexcept(noexcept(_journal.source()))
        {
            return _journal.source();
        }

        constexpr source_type source() && noexcept(noexcept(std::move(_journal).source()))
        {
            return std::move(_journal).source();
        }
        /// @}

        /// @name Iterators
        /// @{
    public:

        constexpr iterator begin() const noexcept
        {
            return iterator{_journal.begin(), _coverages.begin()};
        }

        constexpr iterator end() const noexcept
        {
            return iterator{_journal.end(), _coverages.end()};
        }
        /// @}

        /// @name Modifier
        /// @{
    public:

        template <typename concrete_sequence_t>
            requires std::convertible_to<concrete_sequence_t, sequence_type>
        constexpr iterator record(breakpoint_type breakpoint,
                                  concrete_sequence_t && sequence,
                                  coverage_type coverage)
        {
            // TODO: check if coverage is valid in the sense that it is from the same domain?
            // but how can we ensure this -> by wrapping it and not storing the domain information inside of a coverage object.
            _coverages.reserve(_coverages.size() + 1);
            auto it = _journal.record(std::move(breakpoint), std::forward<concrete_sequence_t>(sequence));
            auto cov_it = _coverages.insert(std::next(_coverages.begin(), std::distance(_journal.begin(), it)),
                                            std::move(coverage));
            return iterator{std::move(it), std::move(cov_it)};
        }
        /// @}

        /// @name Capacities
        /// @{
    public:
        constexpr size_t size() const noexcept
        {
            return _journal.size();
        }

        constexpr size_t max_size() const noexcept
        {
            return _journal.max_size();
        }

        constexpr bool empty() const noexcept
        {
            return _journal.empty();
        }
        /// @}
    };

    template <typename source_t>
    coverage_augmented_breakpoint_multijournal(source_t) -> coverage_augmented_breakpoint_multijournal<source_t>;

    template <libjst::reference_sequence source_t>
    template <typename record_t, typename coverage_t>
    class coverage_augmented_breakpoint_multijournal<source_t>::record_impl
    {
        /// @name Member variables
        /// @{
    private:
        record_t _record;
        coverage_t _coverage;
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr record_impl()
            requires std::default_initializable<record_t> && std::default_initializable<coverage_t> = default;

        constexpr record_impl(record_t record, coverage_t coverage)
            : _record{std::forward<record_t>(record)}, _coverage{std::forward<coverage_t>(coverage)}
        {}

        constexpr sequence_type sequence() const noexcept
        {
            return _record.sequence();
        }

        constexpr coverage_type coverage() const noexcept
        {
            return _coverage;
        }
        /// @}

        /// @name Non-member functions
        /// @{
    public:

        template <typename tag_t, typename self_t>
            requires std::same_as<std::remove_cvref_t<self_t>, record_impl> &&
                     std::invocable<tag_t, record_t>
        constexpr friend auto tag_invoke(tag_t const & tag, self_t && self)
            noexcept(std::is_nothrow_invocable_v<tag_t, record_t>)
            -> std::invoke_result_t<tag_t, record_t>
        {
            return tag(std::forward<record_t>(self._record));
        }

        constexpr friend bool operator==(record_impl const & lhs, record_impl const & rhs) noexcept
        {
            return lhs._record == rhs._record;
        }

        constexpr friend std::weak_ordering operator<=>(record_impl const & lhs, record_impl const & rhs) noexcept
        {
            return lhs._record <=> rhs._record;
        }

        /// @}
    };

    template <libjst::reference_sequence source_t>
    template <bool is_const>
    class coverage_augmented_breakpoint_multijournal<source_t>::iterator_impl
    {
        friend class coverage_augmented_breakpoint_multijournal;

        template <bool>
        friend class iterator_impl;

        /// @name Member types
        /// @{
    private:

        template <typename t>
        using maybe_const_t = std::conditional_t<is_const, t const, t>;

        using journal_iterator = std::ranges::iterator_t<maybe_const_t<base_journal_t>>;
        using coverages_iterator = std::ranges::iterator_t<maybe_const_t<coverages_type>>;
        using base_record_t = std::iter_reference_t<journal_iterator>;
        using coverage_t = std::iter_reference_t<coverages_iterator>;
    public:

        using value_type = record_impl<base_record_t, coverage_t>;
        using reference = value_type;
        using pointer = void;
        using difference_type = std::iter_difference_t<journal_iterator>;
        using iterator_category = std::bidirectional_iterator_tag;
        /// @}

        /// @name Member variables
        /// @{
    private:
        journal_iterator _journal_it{};
        coverages_iterator _coverages_it{};
        /// @}

        /// @name Member functions
        /// @{
    private:

        constexpr explicit iterator_impl(journal_iterator journal_it, coverages_iterator coverages_it)
            noexcept(std::is_nothrow_move_constructible_v<journal_iterator> &&
                     std::is_nothrow_move_constructible_v<coverages_iterator>)
            : _journal_it{std::move(journal_it)}, _coverages_it{std::move(coverages_it)}
        {}
    public:

        constexpr iterator_impl() = default;

        template <bool other_const>
        constexpr iterator_impl(iterator_impl<other_const> other)
            requires (is_const && !other_const)
            : _journal_it{std::move(other._journal_it)}, _coverages_it{std::move(other._coverages_it)}
        {}

        constexpr reference operator*() const noexcept
        {
            return {*_journal_it, *_coverages_it};
        }
        /// @}

        /// @name Arithmetic operators
        /// @{
        constexpr iterator_impl & operator++() noexcept
        {
            ++_journal_it;
            ++_coverages_it;
            return *this;
        }

        constexpr iterator_impl operator++(int) noexcept
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        constexpr iterator_impl & operator--() noexcept
        {
            --_journal_it;
            --_coverages_it;
            return *this;
        }

        constexpr iterator_impl operator--(int) noexcept
        {
            auto tmp = *this;
            --(*this);
            return tmp;
        }
        /// @}

        /// @name Non-member functions
        /// @{

        constexpr friend difference_type operator-(iterator_impl const & lhs, iterator_impl const & rhs) noexcept
        {
            return lhs._journal_it - rhs._journal_it;
        }

        constexpr friend bool operator==(iterator_impl const & lhs, iterator_impl const & rhs) noexcept
        {
            return lhs._journal_it == rhs._journal_it;
        }
        /// @}
    };
}  // namespace libjst
