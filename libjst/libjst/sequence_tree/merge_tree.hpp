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

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>
#include <libjst/variant/breakpoint.hpp>

namespace libjst
{
    template <typename base_tree_t>
    class merge_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;

        class node_impl;

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
            base_node_type base_root = libjst::root(_wrappee);
            auto root_low = base_root.low_boundary();
            return node_impl{std::move(base_root), std::move(root_low)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
   };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::node_impl : public base_node_type {
    public:
        using typename base_node_type::position_type;
    private:

        friend merge_tree_impl;

        position_type _low_boundary{};

        explicit constexpr node_impl(base_node_type && base_node, position_type cached_low) noexcept :
            base_node_type{std::move(base_node)},
            _low_boundary{std::move(cached_low)}
        {}

    public:

        constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit_next<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit_next<false>(base_node_type::next_ref());
        }

        constexpr position_type const & low_boundary() const {
            return _low_boundary;
        }

    private:

        template <bool is_alt_child>
        constexpr std::optional<node_impl> visit_next(auto maybe_child) const {
            if (maybe_child) {
                position_type cached_low = maybe_child->low_boundary();
                node_impl new_child{std::move(*maybe_child), std::move(cached_low)};
                new_child.extend();
                if constexpr (is_alt_child) {
                    new_child.activate_state(node_state::variant);
                }
                return new_child;
            } else {
                return std::nullopt;
            }
        }

        constexpr void extend() {
            while (!base_node_type::high_boundary().is_low_end()) {
                if (auto successor = base_node_type::next_ref(); successor) {
                    static_cast<base_node_type &>(*this) = std::move(*successor);

                } else {
                    break;
                }
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    // template <typename base_tree_t>
    // class merge_tree_impl<base_tree_t>::node_impl::label_impl : public libjst::tree_label_t<base_tree_t> {
    // private:

    //     using base_label_t = libjst::tree_label_t<base_tree_t>;

    //     friend merge_tree_impl;

    //     breakpoint _low_breakend{};
    //     breakpoint _high_breakend{};

    //     template <typename base_label_t>
    //     explicit constexpr label_impl(base_label_t base_label,
    //                                   breakpoint low_breakend,
    //                                   breakpoint high_breakend) noexcept :
    //         base_label_t{std::move(base_label)},
    //         _low_breakend{low_breakend},
    //         _high_breakend{high_breakend}
    //     {}

    // public:

    //     label_impl() = default;

    //     constexpr auto sequence() const noexcept {
    //         assert(_low_breakend <= _high_breakend);
    //         return base_label_t::sequence(_low_breakend.value(), _high_breakend.value());
    //     }
    // };

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
