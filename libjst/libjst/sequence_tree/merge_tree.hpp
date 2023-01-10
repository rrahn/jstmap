// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides merge tree with maximal branch-free nodes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/variant/breakpoint.hpp>

namespace libjst
{
    template <typename base_tree_t>
    class merge_tree_impl {
    private:
        // using wrappee_t = jst::contrib::copyable_box<base_tree_t>;
        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;

        base_tree_t _wrappee{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr merge_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<wrapped_tree_t, merge_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr merge_tree_impl(wrapped_tree_t && wrappee) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee)};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }
   };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend merge_tree_impl;

        class label_impl;

        breakpoint _left_breakpoint{};

        explicit constexpr node_impl(base_node_type && base_node) noexcept :
            base_node_type{std::move(base_node)}
        {
            _left_breakpoint = base_node_type::left_breakpoint();
        }

    public:

        constexpr node_impl() = default;
        constexpr node_impl(node_impl const &) = default;
        constexpr node_impl(node_impl &&) = default;
        constexpr node_impl & operator=(node_impl const &) = default;
        constexpr node_impl & operator=(node_impl &&) = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit(base_node_type::next_ref());
        }

        constexpr auto operator*() const noexcept {
            return label_impl{*static_cast<base_node_type const &>(*this),
                              left_breakpoint(),
                              base_node_type::right_breakpoint()};
        }

    protected:

        constexpr breakpoint left_breakpoint() const noexcept {
            return _left_breakpoint;
        }
    private:

        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                node_impl new_child{std::move(*maybe_child)};
                new_child.extend();
                return new_child;
            } else {
                return std::nullopt;
            }
        }

        // TODO! must be function of base class
        constexpr bool is_nil() const noexcept {
            return base_node_type::from_reference() &&
                //    base_node_type::right_variant() == base_node_type::sink() && // should be the position now?
                   libjst::position(base_node_type::right_variant()) == std::ranges::size(base_node_type::rcs_store().source()) &&
                   base_node_type::get_second_breakpoint_id() != node_descriptor_id::second_first_right;
        }

        constexpr void extend() {
            while (!(is_nil() || base_node_type::is_branching())) {
                if (auto successor = base_node_type::next_ref(); successor) {
                    static_cast<base_node_type &>(*this) = std::move(*successor);
                } else {
                    break;
                }
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::node_impl::label_impl : public libjst::tree_label_t<base_tree_t> {
    private:

        using base_label_t = libjst::tree_label_t<base_tree_t>;

        friend node_impl;

        breakpoint _left_breakpoint{};
        breakpoint _right_breakpoint{};

        template <typename base_label_t>
        explicit constexpr label_impl(base_label_t base_label,
                                      breakpoint left_breakpoint,
                                      breakpoint right_breakpoint) noexcept :
            base_label_t{std::move(base_label)},
            _left_breakpoint{left_breakpoint},
            _right_breakpoint{right_breakpoint}
        {}

    public:

        label_impl() = default;

        constexpr auto sequence() const noexcept {
            assert(_left_breakpoint <= _right_breakpoint);
            return base_label_t::sequence(_left_breakpoint.value(), _right_breakpoint.value());
        }
    };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::sink_impl {
    private:
        friend merge_tree_impl;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    public:
        sink_impl() = default;
    };

    namespace _tree_adaptor {
        inline constexpr struct _merge
        {
            template <typename covered_tree_t, typename ...args_t>
            constexpr auto operator()(covered_tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            merge_tree_impl<std::remove_reference_t<covered_tree_t>>, args_t...>)
                -> merge_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>
            {
                using adapted_tree_t = merge_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>;
                return adapted_tree_t{(covered_tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_merge, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_merge{}, (args_t &&)args...);
            }
        } merge{};
    } // namespace _tree_adaptor

    using _tree_adaptor::merge;
}  // namespace libjst