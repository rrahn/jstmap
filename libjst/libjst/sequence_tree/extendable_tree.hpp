// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>

namespace libjst
{
    template <typename base_tree_t, template <typename, typename> typename node_extension_t>
    class extendable_tree {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;

        using tree_box_t = jst::contrib::copyable_box<base_tree_t>;

        tree_box_t _wrappee{};

    public:

        template <typename wrappee_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, extendable_tree> &&
                      std::constructible_from<tree_box_t, wrappee_t>)
        constexpr explicit extendable_tree(wrappee_t && wrappee) noexcept : _wrappee{(wrappee_t &&)wrappee}
        {}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(*_wrappee)};
        }
        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(*_wrappee)};
        }
    };

    template <typename base_tree_t, template <typename, typename> typename node_extension_t>
    class extendable_tree<base_tree_t, node_extension_t>::node_impl : public base_node_type,
                                                    public node_extension_t<node_impl, base_node_type> {
    private:
        using extension_t = node_extension_t<node_impl, base_node_type>;

        friend extension_t;
        friend extendable_tree;

        explicit constexpr node_impl(base_node_type base_node) noexcept :
            base_node_type{std::move(base_node)}
        {
            extension_t::initialise();
        }

       explicit constexpr node_impl(base_node_type base_node, extension_t extension) noexcept :
            base_node_type{std::move(base_node)},
            extension_t{std::move(extension)}
        {}

    public:

        constexpr std::optional<node_impl> next_alt() const {
            return visit(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const {
            return visit(base_node_type::next_ref());
        }

    private:

        template <typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                extension_t child_extension = extension_t::notify(*maybe_child);
                return node_impl{std::move(*maybe_child), std::move(child_extension)};
            } else {
                return std::nullopt;
            }
        }

        friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t, template <typename, typename> typename node_extension_t>
    class extendable_tree<base_tree_t, node_extension_t>::sink_impl {
    private:
        friend extendable_tree;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

    private:

        sink_impl() = default;

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    };
}  // namespace libjst
