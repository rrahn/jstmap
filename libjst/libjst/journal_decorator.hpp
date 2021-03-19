// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <set>
#include <type_traits>
#include <span>

#include <libjst/journal_entry.hpp>

namespace libjst
{

template <std::ranges::view segment_t>
class journal_decorator
{
private:

    using journal_entry_type = detail::journal_entry<segment_t>;
    using dictionary_type = std::set<journal_entry_type, std::less<void>>;
    using dictionary_iterator = typename dictionary_type::iterator;

    dictionary_type _dictionary{};
    size_t _size{};

public:

    template <bool is_const_range>
    class iterator;

    using value_type = std::ranges::range_value_t<segment_t>;
    using size_type = typename journal_entry_type::size_type;
    using reference = std::ranges::range_reference_t<segment_t>;
    using segment_type = segment_t;

    journal_decorator() = default;
    journal_decorator(journal_decorator const &) = default;
    journal_decorator(journal_decorator &&) = default;
    journal_decorator & operator=(journal_decorator const &) = default;
    journal_decorator & operator=(journal_decorator &&) = default;
    ~journal_decorator() = default;

    explicit journal_decorator(segment_type initial_segment) noexcept
    {
        _size = std::ranges::size(initial_segment);
        if (_size > 0) // only add entry if initial segment is not empty.
            _dictionary.emplace(0, std::move(initial_segment));
    }

    bool record_insertion(size_type const position, segment_type segment)
    {
        if (position > size())
            return false;
        else if (empty() && position != 0)
            return false;

        dictionary_iterator entry_it = _dictionary.lower_bound(position);

        bool inserted = true;
        if (entry_it == _dictionary.end()) // The dictionary is empty and we can just insert the segment.
        {
            inserted = _dictionary.emplace(position, std::move(segment)).second;
        }
        else
        {
            entry_it = split_at(entry_it, position);
            update_remaining_segment_positions(entry_it, std::ranges::size(segment));
            _dictionary.emplace_hint(entry_it, position, std::move(segment));
        }

        update_size(inserted * std::ranges::size(segment));

        assert(check_consistent_segments());
        return inserted;
    }

    bool record_deletion(size_type const first_position, size_type const last_position)
    {
        if (!check_valid_range(first_position, last_position))
            return false;

        int32_t const deletion_size = first_position - last_position;

        update_remaining_segment_positions(erase(first_position, last_position), deletion_size);
        update_size(deletion_size);
        assert(check_consistent_segments());
        return true;
    }

    bool record_substitution(size_type const position, segment_type segment)
    {
        size_type const segment_size = std::ranges::size(segment);
        size_type const last_position = position + segment_size;

        if (!check_valid_range(position, last_position))
            return false;

        _dictionary.emplace_hint(erase(position, last_position), position, segment);

        assert(check_consistent_segments());
        return true;
    }

    size_type size() const noexcept
    {
        return _size;
    }

    bool empty() const noexcept
    {
        return _dictionary.empty();
    }

    iterator<false> begin() noexcept
    {
        return iterator<false>{this, false};
    }

    iterator<true> begin() const noexcept
    {
        return iterator<true>{this, false};
    }

    iterator<false> end() noexcept
    {
        return iterator<false>{this, true};
    }

    iterator<true> end() const noexcept
    {
        return iterator<true>{this, true};
    }

private:

    bool check_valid_range(size_type const first_position, size_type const last_position) const noexcept
    {
        return first_position < last_position && last_position <= size();
    }

    void update_size(int32_t const event_size) noexcept
    {
        _size += event_size;
    }

    auto split_at(std::ranges::iterator_t<dictionary_type> entry_it, size_type const split_position)
    {
        assert(entry_it != _dictionary.end());
        assert(entry_it->segment_begin_position() <= split_position);
        assert(split_position <= entry_it->segment_end_position());

        if (split_position == entry_it->segment_begin_position())
            return entry_it;
        else if (split_position == entry_it->segment_end_position())
            return std::ranges::next(entry_it);

        size_type const right_entry_size = entry_it->segment_end_position() - split_position;
        assert(right_entry_size <= entry_it->segment_size());

        journal_entry_type right_entry{split_position, entry_it->segment().last(right_entry_size)};
        entry_it->segment() = entry_it->segment().first(entry_it->segment_size() - right_entry_size);
        return _dictionary.insert(std::ranges::next(entry_it), std::move(right_entry));
    }

    bool check_consistent_segments() const noexcept
    {
        bool is_consistent{true};
        size_type last_end = 0;
        std::ranges::for_each(_dictionary, [&] (auto const & entry)
        {
            if (entry.segment_begin_position() != last_end)
                is_consistent = false;

            last_end += entry.segment_size();
        });
        return is_consistent;
    }

    void update_remaining_segment_positions(dictionary_iterator first, int32_t const offset)
    {
        std::ranges::for_each(first, _dictionary.end(), [offset] (journal_entry_type const & entry)
        {
            entry.segment_begin_position() += offset;
        });
    }

    dictionary_iterator erase(size_type const first_position, size_type const last_position)
    {
        dictionary_iterator first_entry_it = _dictionary.lower_bound(first_position);

        assert(first_entry_it->segment_begin_position() <= first_position);
        assert(first_position <= first_entry_it->segment_end_position());

        dictionary_iterator first_entry_right_split_it = split_at(first_entry_it, first_position);

        if (first_entry_right_split_it->segment_end_position() >= last_position)
        { // remove infix of current entry.
            size_type const count = last_position - first_position;
            assert(count <= first_entry_right_split_it->segment_size());

            if (count == first_entry_right_split_it->segment_size())
                return _dictionary.erase(first_entry_right_split_it);

            segment_type & right_split_segment = first_entry_right_split_it->segment();
            size_t const remaining_size = first_entry_right_split_it->segment_size() - count;
            first_entry_right_split_it->segment() = right_split_segment.last(remaining_size);
            first_entry_right_split_it->segment_begin_position() += count;
            return first_entry_right_split_it;
        }
        else // delete multiple entries
        {  // remove the suffix of the first entry and the prefix of the last entry.
            dictionary_iterator last_entry_it = _dictionary.lower_bound(last_position);
            dictionary_iterator last_entry_right_split_it = split_at(last_entry_it, last_position);

            assert(*first_entry_right_split_it < *last_entry_right_split_it);

            return _dictionary.erase(first_entry_right_split_it, last_entry_right_split_it);
        }
    }
};

template <std::ranges::contiguous_range range_t>
journal_decorator(range_t &&) -> journal_decorator<std::remove_reference_t<range_t>>;

template <std::ranges::view segment_t>
template <bool is_const_range>
class journal_decorator<segment_t>::iterator
{
private:

    using maybe_const_segment_t = std::conditional_t<is_const_range, segment_t const, segment_t>;

    using segment_iterator = std::ranges::iterator_t<maybe_const_segment_t>;
    using dictionary_iterator = std::ranges::iterator_t<dictionary_type>;

    journal_decorator const * _host{};
    segment_iterator _segment_it{};
    dictionary_iterator _entry_it{};

    template <bool>
    friend class iterator;

public:

    using value_type = std::ranges::range_value_t<maybe_const_segment_t>;
    using reference = std::ranges::range_reference_t<maybe_const_segment_t>;
    using pointer = std::add_pointer_t<value_type>;
    using difference_type = std::ranges::range_difference_t<maybe_const_segment_t>;
    using iterator_category = std::random_access_iterator_tag; // we should not but well.

    iterator() = default;
    iterator(iterator const &) = default;
    iterator(iterator &&) = default;
    iterator & operator=(iterator const &) = default;
    iterator & operator=(iterator &&) = default;
    ~iterator() = default;

    iterator(journal_decorator const * host, bool const is_end) : _host{host}
    {
        _entry_it = (is_end) ? _host->_dictionary.end() : _host->_dictionary.begin();

        if (!_host->_dictionary.empty())
            _segment_it = (!end_of_journal()) ? segment_begin() : set_end();
    }

    iterator(iterator<!is_const_range> other)
        requires is_const_range
    :
        _host{other._host},
        _segment_it{other._segment_it},
        _entry_it{other._entry_it}
    {}

    /*!\name Element access
     * \{
     */
    reference operator*() const noexcept
    {
        return *_segment_it;
    }

    pointer operator->() const noexcept
    {
        return std::addressof(*_segment_it);
    }

    reference operator[](difference_type const offset) const noexcept
    {
        return *(*this + offset);
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    iterator & operator++() noexcept
    {
        if (++_segment_it == segment_end())
        {
            ++_entry_it;
            _segment_it = (!end_of_journal()) ? segment_begin() : set_end();
        }

        return *this;
    }

    iterator operator++(int) noexcept(std::is_nothrow_copy_constructible_v<iterator>)
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    iterator & operator+=(difference_type offset) noexcept
    {
        // if offset is negative, use the preceding elements of current segment to check if binary search is needed.
        // if offset is positive use remaining elements of current segment to check if binary search is needed.
        segment_iterator segment_terminus = (offset < 0) ? segment_begin() : segment_end();
        if (std::abs(std::ranges::distance(_segment_it, segment_terminus)) <= std::abs(offset))
        {
            difference_type next_position = segment_position() + offset;
            _entry_it = _host->_dictionary.upper_bound(next_position);

            assert([&] { return end_of_journal() || next_position >= _entry_it->segment_begin_position(); }());
            _segment_it = (!end_of_journal())
                        ? std::ranges::next(segment_begin(), next_position - _entry_it->segment_begin_position())
                        : set_end();
        }
        else
        {
            std::ranges::advance(_segment_it, offset);
        }

        return *this;
    }

    iterator operator+(difference_type const offset) const noexcept
    {
        iterator tmp{*this};
        return tmp += offset;
    }

    friend iterator operator+(difference_type const offset, iterator const & rhs) noexcept
    {
        return rhs + offset;
    }

    iterator & operator--() noexcept
    {
        if (_segment_it == segment_begin())
        {
            int32_t offset = -(_entry_it != _host->_dictionary.begin());
            std::advance(_entry_it, offset);
            _segment_it = (offset == 0) ? _segment_it : std::ranges::prev(segment_end());
        }
        else
        {
            --_segment_it;
        }
        return *this;
    }

    iterator operator--(int) noexcept(std::is_nothrow_copy_constructible_v<iterator>)
    {
        iterator tmp{*this};
        --(*this);
        return tmp;
    }

    iterator & operator-=(difference_type offset) noexcept
    {
        return (*this) += -offset;
    }

    iterator operator-(difference_type offset) const noexcept
    {
        iterator tmp{*this};
        return tmp -= offset;
    }

    template <bool other_const_range>
    difference_type operator-(iterator<other_const_range> const & rhs) const noexcept
    {
        return segment_position() - rhs.segment_position();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    template <bool other_const_range>
    bool operator==(iterator<other_const_range> const & rhs) const noexcept
    {
        return (_entry_it == rhs._entry_it) && (_segment_it == rhs._segment_it);
    }

    template <bool other_const_range>
    std::strong_ordering operator<=>(iterator<other_const_range> const & rhs) const noexcept
    {
        return segment_position() <=> rhs.segment_position();
    }
    //!\}

private:

    bool end_of_journal() const noexcept
    {
        return _entry_it == _host->_dictionary.end();
    }

    auto set_end() noexcept
    {
        assert(end_of_journal());
        assert(_entry_it != _host->_dictionary.begin());

        --_entry_it;
        return segment_end();
    }

    auto segment_begin() const noexcept
    {
        assert(_host != nullptr);
        assert(_entry_it != _host->_dictionary.end());

        return std::ranges::begin(_entry_it->segment());
    }

    auto segment_end() const noexcept
    {
        assert(_host != nullptr);
        assert(_entry_it != _host->_dictionary.end());

        return std::ranges::end(_entry_it->segment());
    }

    difference_type segment_position() const noexcept
    {
        assert(_host != nullptr);

        return (_entry_it != _host->_dictionary.end())
            ? _entry_it->segment_begin_position() + std::ranges::distance(segment_begin(), _segment_it)
            : _host->size();
    }

};

}  // namespace libjst
