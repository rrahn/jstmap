// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cstring>
#include <utility>
#include <vector>
#include <tuple>

#include <seqan/vcf_io.h>

#include <libjst/coverage/range_domain.hpp>

#include <jstmap/global/jstmap_types.hpp>

namespace jstmap
{
    struct variant_stat{
        size_t snv_count{};
        size_t indel_count{};
    };
    class stripped_vcf_record
    {
    private:
        using alternative_t = std::string;
        using genotypes_t = std::vector<coverage_t>;
        using position_t = uint32_t;
        using coverage_value_t = typename coverage_t::value_type;
        using domain_t = libjst::range_domain<coverage_value_t>;

        // relevant fields to set.
        std::string _ref{};
        std::string _chrom_name{};
        std::vector<alternative_t> _alt{};
        genotypes_t _genotypes{};
        domain_t _domain{};
        variant_stat * _stat{};
        position_t _pos{};
        int32_t _sample_count{};
        int32_t _haplotype_count{};
        int16_t _chrom_id{};
        uint8_t _alternative_count{};

    public:

        using intermediate_variant_t = std::tuple<int32_t, std::string, int32_t, coverage_t>;

        stripped_vcf_record() = default;

        template <typename vcf_file_t>
        stripped_vcf_record(vcf_file_t & vcf_file, domain_t domain, variant_stat & stat) :
            _domain{std::move(domain)},
            _stat{&stat}
        {
            auto & file_context = seqan2::context(vcf_file);
            _sample_count = seqan2::length(seqan2::sampleNames(file_context));
            _haplotype_count = _sample_count << 1; // TODO: detect ploidy
            read_record(file_context, directionIterator(vcf_file, seqan2::Input{}));
            _chrom_id = nameToId(contigNamesCache(file_context), seqan2::CharString{_chrom_name});
        }

        template <typename vcf_context_t>
        std::string const & contig_name(vcf_context_t const & vcf_context) const noexcept
        {
            return _chrom_name;
        }

        genotypes_t const & field_genotype() const noexcept;
        void alternatives(rcs_store_t &) const;

    private:
        void set_field_chrom(std::string_view);

        void set_field_pos(std::string_view);

        void set_field_ref(std::string_view);

        void set_field_alt(std::string_view);

        void set_field_genotype(std::string_view);

        std::string_view read_field(std::string_view &) noexcept;

        template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
        inline void
        read_record(seqan2::VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context, TForwardIter & iter)
        {
            if (_sample_count == 0)
                return;

            // get the next line on the buffer.
            clear(context.buffer);
            readLine(context.buffer, iter);
            std::string_view buffer{toCString(context.buffer), toCString(context.buffer) + length(context.buffer)};
            // Parse field #CHROM
            set_field_chrom(read_field(buffer));
            // Parse field #POS
            set_field_pos(read_field(buffer));
            // Skip field #ID -- to annotate variants, e.g. dbSNP identifier
            read_field(buffer);
            // Parse field #REF
            set_field_ref(read_field(buffer));
            // Parse field #ALT
            set_field_alt(read_field(buffer));
            // Skip field #QUAL
            read_field(buffer);
            // Skip field #FILTER
            read_field(buffer);
            // Skip field #INFO
            read_field(buffer);
            // Skip field #FORMAT
            read_field(buffer);
            // Parse genotypes
            set_field_genotype(buffer);
        }
    };
}  // namespace jstmap
