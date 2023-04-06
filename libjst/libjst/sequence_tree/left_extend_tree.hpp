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
#include <libjst/sequence_tree/breakend_site_trimmed.hpp>

namespace libjst
{

    template <typename base_tree_t>
    class left_extend_tree_impl {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using offset_type = std::ptrdiff_t;

        class node_impl;

        base_tree_t _wrappee{};
        offset_type _offset{};

    public:

        template <typename wrappee_t, std::integral offset_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, left_extend_tree_impl> &&
                      std::constructible_from<base_tree_t, wrappee_t>)
        constexpr explicit left_extend_tree_impl(wrappee_t && wrappee, offset_t offset) noexcept :
            _wrappee{(wrappee_t &&)wrappee},
            _offset{static_cast<offset_type>(offset)}
        {
            assert(_offset >= 0);
        }

        constexpr node_impl root() const noexcept {
            base_node_type base_root = libjst::root(_wrappee);
            offset_type lowest = libjst::position(base_root.low_boundary());
            return node_impl{libjst::root(_wrappee), _offset, lowest};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
    };

    template <typename base_tree_t>
    class left_extend_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:
        using base_position_type = typename base_node_type::position_type;

        friend left_extend_tree_impl;

        offset_type _offset{};
        offset_type _lowest{};

        explicit constexpr node_impl(base_node_type && base_node, offset_type offset, offset_type lowest) noexcept :
            base_node_type{std::move(base_node)},
            _offset{offset},
            _lowest{lowest}
        {}

    public:

        using position_type = breakend_site_trimmed<base_position_type>;

        node_impl() = default;

        // constexpr auto operator*() const noexcept {
        //     return label_impl{*static_cast<base_node_type const &>(*this), _offset, _lowest};
        // }

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit(base_node_type::next_ref());
        }

        constexpr position_type low_boundary() const {
            using position_value_t = typename position_type::position_value_type;
            base_position_type const & base_boundary = base_node_type::low_boundary();
            position_value_t low_position = std::max<offset_type>(libjst::position(base_boundary) - _offset, _lowest);
            return position_type{base_boundary, low_position};
        }

        constexpr position_type high_boundary() const {
            return position_type{base_node_type::high_boundary()};
        }

    private:

        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                return node_impl{std::move(*maybe_child), _offset, _lowest};
            } else {
                return std::nullopt;
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    // template <typename base_tree_t>
    // class left_extend_tree_impl<base_tree_t>::label_impl : public base_label_type {
    // private:

    //     friend left_extend_tree_impl;

    //     offset_type _left_extension{};
    //     offset_type _min_position{};

    //     constexpr explicit label_impl(base_label_type base_label, offset_type left_extension, offset_type min_position) noexcept :
    //         base_label_type{std::move(base_label)},
    //         _left_extension{left_extension},
    //         _min_position{min_position}
    //     {}

    // public:

    //     using typename base_label_type::size_type;

    //     label_impl() = default;

    //     constexpr auto sequence(size_type left_pos, size_type right_pos = base_label_type::npos) const {
    //         assert(left_pos >= _min_position);
    //         assert(left_pos <= right_pos);
    //         left_pos = std::max<std::ptrdiff_t>(_min_position, left_pos - _left_extension);
    //         return base_label_type::sequence(left_pos, right_pos);
    //     }

    //     constexpr auto sequence() const {
    //         return sequence(_min_position);
    //     }
    // };

    namespace _tree_adaptor {
        inline constexpr struct _left_extend
        {
            template <typename labelled_tree_t, std::integral left_extension_t>
            constexpr auto operator()(labelled_tree_t && tree, left_extension_t left_extension) const
                -> left_extend_tree_impl<labelled_tree_t>
            {
                assert(left_extension >= 0);
                return left_extend_tree_impl<labelled_tree_t>{(labelled_tree_t &&) tree, std::move(left_extension)};
            }

            template <std::integral left_extension_t>
            constexpr auto operator()(left_extension_t const left_extension) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, left_extension_t>)
                -> jst::contrib::closure_result_t<_left_extend, left_extension_t>
            { // we need to store the type that needs to be called later!
                assert(left_extension >= 0);
                return jst::contrib::make_closure(_left_extend{}, left_extension);
            }
        } left_extend{};
    } // namespace _tree_adaptor

    using _tree_adaptor::left_extend;
}  // namespace libjst
