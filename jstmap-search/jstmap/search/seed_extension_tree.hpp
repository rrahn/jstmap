// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides extension base tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libspm/closure_adaptor.hpp>
#include <libspm/std/tag_invoke.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/breakend_site_min.hpp>
#include <libjst/sequence_tree/breakend_site_trimmed.hpp>
#include <libjst/variant/concept.hpp>

#include <jstmap/global/match_position.hpp>

namespace jstmap
{
    template <typename wrapped_tree_t>
    class seed_extension_tree
    {
    private:
        using base_node_type = libjst::tree_node_t<wrapped_tree_t>;
        using sink_type = libjst::tree_sink_t<wrapped_tree_t>;
        using base_cargo_type = libjst::tree_label_t<wrapped_tree_t>;
        using difference_type = std::ptrdiff_t;

        // generate local position!

        class node_impl;
        class cargo_impl;

        wrapped_tree_t _wrappee{};
        base_node_type _base{};
        difference_type _min_low_position{};
        difference_type _max_label_size{};

    public:

        template <typename wrappee_t, std::unsigned_integral extension_size_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, seed_extension_tree>) &&
                      std::constructible_from<wrapped_tree_t, wrappee_t>
        seed_extension_tree(wrappee_t && wrappee, match_position start, extension_size_t max_label_size) noexcept :
            _wrappee{(wrappee_t &&) wrappee}
        {
            _base = _wrappee.seek(start.tree_position); // from seeking to the node!
            auto base_cargo = *_base;
            // - compute the label distance until next breakend.
            // log_debug("std::ranges::ssize(base_cargo.path_sequence()): ", std::ranges::ssize(base_cargo.path_sequence()));
            // log_debug("position(*(_base.low_boundary())): ", libjst::position(*(_base.low_boundary())));
            // log_debug("position(*(_base.high_boundary())): ", libjst::position(*(_base.high_boundary())));
            assert(std::ranges::ssize(base_cargo.path_sequence()) >= start.label_offset);
            difference_type _distance_to_high = std::ranges::ssize(base_cargo.path_sequence()) - start.label_offset;

            _min_low_position = libjst::position(_base.high_boundary()) - _distance_to_high;

            assert(static_cast<difference_type>(libjst::position(_base.low_boundary())) <= _min_low_position);
            // - update remaining max label size of root node.
            _max_label_size = static_cast<difference_type>(max_label_size) - _distance_to_high;
        }

        constexpr node_impl root() const noexcept {
            return node_impl{_base, _min_low_position, _max_label_size};
        }

        constexpr sink_type sink() const noexcept {
            return _wrappee.sink();
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
    };

    template <typename wrapped_tree_t>
    class seed_extension_tree<wrapped_tree_t>::node_impl : public base_node_type
    {
    private:
        friend seed_extension_tree;

        using base_low_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().low_boundary())>;
        using base_high_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().high_boundary())>;
        using position_value_type = libjst::variant_position_t<base_low_position_type>;

        difference_type _min_low_position{};
        difference_type _remaining_label_size{};

        explicit constexpr node_impl(base_node_type base_node,
                                     difference_type min_low_position,
                                     difference_type remaining_label_size) noexcept :
            base_node_type{std::move(base_node)},
            _min_low_position{std::move(min_low_position)},
            _remaining_label_size{remaining_label_size}
        {}

    public:

        using low_position_type = libjst::breakend_site_min<base_low_position_type>;
        using high_position_type = libjst::breakend_site_trimmed<base_high_position_type>;

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (is_leaf())
                return std::nullopt;
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (is_leaf())
                return std::nullopt;
            return visit<false>(base_node_type::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

        constexpr low_position_type low_boundary() const noexcept {
            return low_position_type{base_node_type::low_boundary(),
                                     static_cast<position_value_type>(_min_low_position)};
        }

        constexpr high_position_type high_boundary() const noexcept {
            auto high_position = base_node_type::high_boundary();
            position_value_type max_position = libjst::position(high_position) + _remaining_label_size;
            return high_position_type{std::move(high_position), max_position};
        }

    private:

        constexpr bool is_leaf() const noexcept {
            return _remaining_label_size <= 0;
        }

        template <bool is_alt_node, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                difference_type next_span = libjst::position(maybe_child->high_boundary()) -
                                            libjst::position(maybe_child->low_boundary());
                if constexpr (is_alt_node) {
                    auto delta = *(maybe_child->low_boundary());
                    next_span += libjst::effective_size(delta) - std::ranges::ssize(libjst::alt_sequence(delta));
                }
                return node_impl{std::move(*maybe_child), _min_low_position, _remaining_label_size - next_span};
            } else {
                return std::nullopt;
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename wrapped_tree_t>
    class seed_extension_tree<wrapped_tree_t>::cargo_impl : public base_cargo_type
    {
    private:
        friend seed_extension_tree;

        node_impl const * _node{};
        constexpr explicit cargo_impl(node_impl const * node) noexcept :
            base_cargo_type{*static_cast<base_node_type const &>(*node)},
            _node{node}
        {
        }

    public:
        constexpr auto sequence() const noexcept {
            assert(_node != nullptr);

            return base_cargo_type::sequence(libjst::position(_node->low_boundary()),
                                             libjst::position(_node->high_boundary()));
        }

        constexpr auto path_sequence() const noexcept {
            assert(_node != nullptr);
            return base_cargo_type::sequence(0, libjst::position(_node->high_boundary()));
        }

        constexpr bool is_leaf() const noexcept {
            assert(_node != nullptr);
            return _node->is_leaf();
        }

    protected:
        using base_cargo_type::sequence;
    };

    namespace _tree_adaptor {
        struct _extend_from
        {
            template <typename tree_t, std::unsigned_integral extension_size_t>
            constexpr auto operator()(tree_t && tree, match_position start_position, extension_size_t extension_size) const
                noexcept(std::is_nothrow_constructible_v<seed_extension_tree<std::remove_reference_t<tree_t>>>)
                -> seed_extension_tree<std::remove_reference_t<tree_t>>
            {
                using adapted_tree_t = seed_extension_tree<std::remove_reference_t<tree_t>>;
                return adapted_tree_t{(tree_t &&)tree, std::move(start_position), extension_size};
            }

            template <std::unsigned_integral extension_size_t>
            constexpr auto operator()(match_position start_position, extension_size_t extension_size) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<spm::make_closure>, _extend_from, match_position, extension_size_t>)
                -> spm::closure_result_t<_extend_from, match_position, extension_size_t>
            {
                return spm::make_closure(_extend_from{}, std::move(start_position), extension_size);
            }
        };
        inline constexpr _extend_from extend_from{};
    } // namespace _tree_adaptor

    using _tree_adaptor::extend_from;
}  // namespace jstmap
