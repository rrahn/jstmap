// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides int coverage.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <libjst/coverage/concept.hpp>
#include <libjst/coverage/range_domain.hpp>
#include <libjst/utility/sorted_vector.hpp>

namespace libjst
{

    template <std::unsigned_integral value_t>
    class int_coverage {

        using coverage_domain_t = range_domain<value_t>;
        using domain_value_type = typename coverage_domain_t::value_type;
        using data_type = sorted_vector<domain_value_type>;

        data_type _data{};
        [[no_unique_address]] coverage_domain_t _domain{};

        // needs its own iterator and output iterator?

    public:
        using value_type = domain_value_type;
        using iterator = std::ranges::iterator_t<data_type const>;

        explicit constexpr int_coverage() = default;
        explicit constexpr int_coverage(coverage_domain_t domain) noexcept :
            _data{},
            _domain{std::move(domain)}
        {}

        template <typename elem_range_t>
            requires (!std::same_as<std::remove_cvref_t<elem_range_t>, int_coverage>) &&
                     (!std::same_as<std::remove_cvref_t<elem_range_t>, std::initializer_list<value_type>>) &&
                      std::integral<std::ranges::range_value_t<elem_range_t>>
        explicit constexpr int_coverage(elem_range_t && from_list, coverage_domain_t domain) :
            int_coverage{std::move(domain)}
        {
            std::ranges::for_each(from_list, [&] (value_type elem) {
                if (!get_domain().is_member(elem))
                    throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};

                _data.emplace_hint(_data.end(), std::move(elem));
            });
        }

        explicit constexpr int_coverage(std::initializer_list<value_type> from_list, coverage_domain_t domain) :
            int_coverage{std::move(domain)}
        {
            std::ranges::for_each(from_list, [&] (value_type elem) {
                if (!get_domain().is_member(elem))
                    throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};

                _data.emplace_hint(_data.end(), std::move(elem));
            });
        }

        constexpr bool contains(std::ptrdiff_t idx) const noexcept {
            return _data.contains(idx);
        }

        constexpr iterator insert(value_type elem) {
            if (!get_domain().is_member(elem))
                throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};

            return _data.insert(elem);
        }

        constexpr iterator insert(iterator hint, value_type elem) {
            if (!get_domain().is_member(elem))
                throw std::domain_error{"The given element " + std::to_string(elem) + " is no member of the coverage domain!"};

            return _data.insert(std::move(hint), elem);
        }

        constexpr void clear() noexcept {
            _data.clear();
        }

        constexpr iterator erase(iterator first) noexcept {
            return erase(first, std::ranges::next(first));
        }

        constexpr iterator erase(iterator first, iterator last) noexcept {
            return _data.erase(first, last);
        }

        constexpr value_type front() const noexcept {
            assert(!empty());
            return *begin();
        }

        constexpr value_type back() const noexcept {
            assert(!empty());
            return *std::ranges::prev(end());
        }

        constexpr bool empty() const noexcept {
            return _data.empty();
        }

        constexpr size_t size() const noexcept {
            return _data.size();
        }

        constexpr size_t max_size() const noexcept {
            return get_domain().size();
        }

        constexpr void reserve(size_t const new_capacity) noexcept {
            return _data.reserve(new_capacity);
        }

        constexpr coverage_domain_t const & get_domain() const noexcept {
            return _domain;
        }

        constexpr bool any() const noexcept {
            return !_data.empty();
        }

        constexpr iterator begin() const noexcept {
            return _data.begin();
        }

        constexpr iterator end() const noexcept {
            return _data.end();
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <typename archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_data, _domain);
        }

        template <typename archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_data, _domain);
        }

    private:

        constexpr friend bool operator==(int_coverage const &, int_coverage const &) noexcept = default;

        constexpr friend int_coverage
        tag_invoke(libjst::tag_t<coverage_intersection>, int_coverage const & first, int_coverage const & second) {
            // if (first.get_domain() != second.get_domain())
            //     throw std::domain_error{"Trying to intersect elements from different coverage domains."};
            // we get the first as copy.
            // we use two iterators to compare the elements with each other.
            // if first < second: move until first >= second
            if (first.size() == first.max_size()) return second;
            if (second.size() == second.max_size()) return first;

            // return first.compute_intersection(second);

            size_t min_size = std::min(first.size(), second.size());
            int_coverage result{};
            result.reserve(min_size);
            std::ranges::set_intersection(first._data.data(),
                                          second._data.data(),
                                          std::back_inserter(result._data.data()));
            return result;
        }

        constexpr friend int_coverage
        tag_invoke(libjst::tag_t<coverage_difference>, int_coverage const & first, int_coverage const & second) {
            // if (first.get_domain() != second.get_domain())
            //     throw std::domain_error{"Trying to intersect elements from different coverage domains."};

            int_coverage result{};
            std::ranges::set_difference(first, second, std::inserter(result, result.end()));
            return result;
        }

        constexpr int_coverage compute_intersection(int_coverage rhs) const noexcept {
            auto lhs_it = _data.data().begin();
            auto rhs_it = rhs._data.data().begin();
            auto inserter = rhs_it;

            while (lhs_it != _data.data().end() && rhs_it != rhs._data.data().end()) {
                if (*lhs_it == *rhs_it) {
                    std::iter_swap(inserter, rhs_it);
                    ++inserter;
                    ++lhs_it;
                    ++rhs_it;
                } else {
                    bool const left_less = *lhs_it < *rhs_it;
                    lhs_it += left_less;
                    rhs_it += !left_less;
                }
            }
            rhs._data.data().erase(inserter, rhs._data.data().end());
            rhs._data.data().shrink_to_fit();
            return rhs;
        }
    };
}  // namespace libjst
