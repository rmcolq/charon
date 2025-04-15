#ifndef CHARON_INDEX_H
#define CHARON_INDEX_H

#pragma once

#include <unordered_map>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <plog/Log.h>

#include <index_main.hpp>
#include <input_summary.hpp>
#include <input_stats.hpp>


class Index
{
    private:
        uint8_t window_size_{};
        uint8_t kmer_size_{};
        double max_fpr_{};
        InputSummary summary_{};
        InputStats stats_{};
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_{};

    public:
        static constexpr uint32_t version{3u};

        Index() = default;
        Index(Index const &) = default;
        Index(Index &&) = default;
        Index & operator=(Index const &) = default;
        Index & operator=(Index &&) = default;
        ~Index() = default;

        Index(const IndexArguments & arguments, const InputSummary & summary, const InputStats & stats, const seqan3::interleaved_bloom_filter<seqan3::uncompressed>& ibf):
                window_size_{arguments.window_size},
                kmer_size_{arguments.kmer_size},
                max_fpr_{arguments.max_fpr},
                summary_{summary},
                stats_{stats},
                ibf_(ibf)
        {}

        uint8_t window_size() const
        {
            return window_size_;
        }

        uint8_t kmer_size() const
        {
            return kmer_size_;
        }

        uint8_t num_bins() const
        {
            return summary_.num_bins;
        }

        uint8_t num_categories() const
        {
            return summary_.num_categories();
        }

        std::vector<std::string> categories() const
        {
            return summary_.categories;
        }

        uint8_t get_host_index() const
        {
            const auto index1 = summary_.category_index("host");
            const auto index2 = summary_.category_index("human");
            auto index = std::min(index1, index2);
            if (index == std::numeric_limits<uint8_t>::max())
                PLOG_ERROR << "Index does not contain 'host' or 'human' as a category ";
            assert(index < std::numeric_limits<uint8_t>::max());
            return index;
        }

        double max_fpr() const
        {
            return max_fpr_;
        }

        InputStats stats() const
        {
            return stats_;
        }

        InputSummary summary() const
        {
            return summary_;
        }

        std::unordered_map<uint8_t, std::string> bin_to_category() const
        {
            return summary_.bin_to_category;
        }

        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> const & ibf() const
        {
            return ibf_;
        }

        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>::membership_agent_type agent() const
        {
            return ibf_.membership_agent();
        }

        /*!\cond DEV
         * \brief Serialisation support function.
         * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
         * \param[in] archive The archive being serialised from/to.
         *
         * \attention These functions are never called directly.
         * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
         */
        template <seqan3::cereal_archive archive_t>
        void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
        {
            try
                {
                    archive(window_size_);
                    archive(kmer_size_);
                    archive(max_fpr_);
                    archive(summary_);
                    archive(stats_);
                    archive(ibf_);
                }
                    // GCOVR_EXCL_START
                catch (std::exception const & e)
                {
                    PLOG_ERROR << "Cannot read index: " + std::string{e.what()};
                    exit(1);
                }

        }

        /* \brief Serialisation support function. Do not load the actual data.
         * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_input_archive.
         * \param[in] archive The archive being serialised from/to.
         *
         * \attention These functions are never called directly.
         * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
         */
        template <seqan3::cereal_input_archive archive_t>
        void load_parameters(archive_t & archive)
        {
            try
                {
                    archive(window_size_);
                    archive(kmer_size_);
                    archive(max_fpr_);
                    archive(summary_);
                    archive(stats_);
                    archive(ibf_);
                }
                    // GCOVR_EXCL_START
                catch (std::exception const & e)
                {
                    PLOG_ERROR << "Cannot read index: " + std::string{e.what()};
                    exit(1);
                }
        }
        //!\endcond
    };

#endif // CHARON_INDEX_H
