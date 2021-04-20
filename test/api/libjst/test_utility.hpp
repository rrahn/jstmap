// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <map>
#include <string_view>
#include <vector>

#include <cereal/archives/binary.hpp>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>

#include <libjst/context_position.hpp>

namespace libjst::test
{

template <typename jst_t>
jst_t load_jst(std::filesystem::path const & jst_file)
{
    using namespace std::literals;

    std::fstream jst_input_stream{jst_file};

    if (!jst_input_stream.good())
        throw std::runtime_error{"Couldn't open path for loading the jst! The path is ["s + jst_file.string() + "]"s};

    jst_t jst{};

    {
        cereal::BinaryInputArchive input_archive{jst_input_stream};
        jst.load(input_archive);
    }

    return jst;
}

inline constexpr auto make_gapped = [] (std::string_view const seq) -> std::vector<seqan3::gapped<char>>
{
    std::vector<seqan3::gapped<char>> tmp{};
    tmp.reserve(seq.size());

    std::for_each(seq.begin(), seq.end(), [&] (char const c)
    {
        if (c == '-')
            tmp.emplace_back(seqan3::gap{});
        else
            tmp.emplace_back(c);
    });

    return tmp;
};

inline constexpr auto sequence_to_string = [] <seqan3::sequence range_t>(range_t && sequence) -> std::string
{
    return sequence | seqan3::views::to_char | seqan3::views::to<std::string>;
};

class jst_context_map_fixture
{
public:
    using context_position_map_t = std::map<std::string_view, std::vector<libjst::context_position>>;

    context_position_map_t context_position_map{};

    // Variable to validate a correct traversal.
    int64_t total_context_count{};
    std::vector<libjst::context_position> unknown_locations{};

    bool all_contexts_enumerated() const
    {
        return total_context_count == 0;
    }

    template <typename position_range_t>
        requires std::same_as<std::ranges::range_value_t<position_range_t>, libjst::context_position>
    bool context_positions_exist(std::string_view context, position_range_t && locations)
    {
        if (std::ranges::empty(locations))
            return true;

        if (auto it = context_position_map.find(context); it != context_position_map.end())
        {
            bool found_all{true};
            for (libjst::context_position const & actual_location : locations)
            {
                size_t erased_elements = std::erase(it->second, actual_location);

                EXPECT_LE(erased_elements, 1u);

                if (erased_elements == 0u)
                {
                    unknown_locations.push_back(actual_location);
                    found_all = false;
                }

                --total_context_count;
            }
            return found_all;
        }
        return  false;
    }

    void print_unvisited_contexts() const
    {
        for (auto && [context, positions] : context_position_map)
        {
            if (positions.empty())
                continue;

            std::cout << "Context: " << context;
            for (auto && [id, pos] : positions)
                std::cout << "\t [" << id << ", " << pos << "]";

            std::cout << "\n";
        }
    }

    void print_unknown_context_locations() const
    {
        for (libjst::context_position const & unkown_location : unknown_locations)
            std::cout << unkown_location << "\n";
    }

protected:

    void generate_context_map(size_t const context_size, std::vector<std::string> const & sequences)
    {
        size_t sequence_index = 0;
        std::ranges::for_each(sequences, [&] (std::string const & sequence)
        {
            std::string_view sv{sequence};
            size_t context_end_position = std::max<int32_t>(sv.size() - context_size + 1, 0);
            assert(context_end_position <= sv.size());
            for (size_t context_position = 0; context_position < context_end_position; ++context_position)
            {
                libjst::context_position context_location{.sequence_id = sequence_index,
                                                          .sequence_position = context_position};
                using value_t = context_position_map_t::value_type;
                value_t insert_value{sv.substr(context_position, context_size), std::vector{context_location}};
                if (auto [it, inserted] = context_position_map.insert(std::move(insert_value)); !inserted)
                    it->second.emplace_back(context_location);

                ++total_context_count;
            }
            ++sequence_index;
        });
    }
};

} // namespace libjst::test
