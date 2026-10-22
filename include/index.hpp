#ifndef CHARON_INDEX_H
#define CHARON_INDEX_H

#pragma once

#include <array>
#include <string>
#include <cinttypes>
#include <optional>
#include <variant>
#include <ankerl/unordered_dense.h>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <plog/Log.h>

#include <index_main.hpp>
#include <input_summary.hpp>
#include <input_stats.hpp>
#include <layout.hpp>


class Index {
public:
    using ibf_uncompressed = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>;
    using ibf_compressed = seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>;

private:
    uint8_t window_size_{};
    uint8_t kmer_size_{};
    double max_fpr_{};
    InputSummary summary_{};
    InputStats stats_{};
    Layout layout_{};
    // Indexes written by Charon <= v3 store a compressed IBF and are queried compressed (as before).
    // Newer indexes store and query the uncompressed IBF, which is much faster to query.
    std::variant<ibf_uncompressed, ibf_compressed> ibf_{};

public:
    static constexpr uint32_t version{4u};
    // Written before the cereal archive in version >= 4 files; older files start straight with the archive.
    static constexpr std::array<char, 8> format_magic{'C', 'H', 'A', 'R', 'O', 'N', 'v', '4'};

    Index() = default;

    Index(Index const &) = default;

    Index(Index &&) = default;

    Index &operator=(Index const &) = default;

    Index &operator=(Index &&) = default;

    ~Index() = default;

    Index(const IndexArguments &arguments, const InputSummary &summary, const InputStats &stats,
          ibf_uncompressed ibf, std::optional<Layout> layout = std::nullopt) :
            window_size_{arguments.window_size},
            kmer_size_{arguments.kmer_size},
            max_fpr_{arguments.max_fpr},
            summary_{summary},
            stats_{stats},
            layout_{layout.value_or(Layout::identity(summary.num_bins))},
            ibf_(std::move(ibf)) {
        layout_.finalise();
    }

    [[nodiscard]] constexpr uint8_t window_size() const {
        return window_size_;
    }

    [[nodiscard]] constexpr uint8_t kmer_size() const {
        return kmer_size_;
    }

    [[nodiscard]] constexpr uint8_t num_bins() const {
        return summary_.num_bins;
    }

    [[nodiscard]] uint8_t num_categories() const {
        return summary_.num_categories();
    }

    [[nodiscard]] const std::vector<std::string>& categories() const {
        return summary_.categories;
    }

    [[nodiscard]] uint8_t get_host_index() const {
        const auto index1 = summary_.category_index("host");
        const auto index2 = summary_.category_index("human");
        auto index = std::min(index1, index2);
        if (index == std::numeric_limits<uint8_t>::max()) {
            PLOG_ERROR << "Index does not contain 'host' or 'human' as a category ";
        }
        assert(index < std::numeric_limits<uint8_t>::max());
        return index;
    }

    [[nodiscard]] uint8_t get_category_index(const std::string category) const {
        const auto index = summary_.category_index(category);
        if (index == std::numeric_limits<uint8_t>::max()) {
            PLOG_ERROR << "Index does not contain category ";
        }
        assert(index < std::numeric_limits<uint8_t>::max());
        return index;
    }

    [[nodiscard]] double max_fpr() const {
        return max_fpr_;
    }

    [[nodiscard]] const InputStats& stats() const {
        return stats_;
    }

    [[nodiscard]] const InputSummary& summary() const {
        return summary_;
    }

    [[nodiscard]] ankerl::unordered_dense::map<uint8_t, std::string> bin_to_category() const {
        return summary_.bin_to_category;
    }

    [[nodiscard]] const Layout &layout() const {
        return layout_;
    }

    [[nodiscard]] bool is_compressed() const {
        return std::holds_alternative<ibf_compressed>(ibf_);
    }

    [[nodiscard]] size_t technical_bins() const {
        return std::visit([](auto const &ibf) { return ibf.bin_count(); }, ibf_);
    }

    [[nodiscard]] size_t bin_size() const {
        return std::visit([](auto const &ibf) { return ibf.bin_size(); }, ibf_);
    }

    // Calls f(agent) with a membership agent for whichever IBF layout this index holds.
    // Call once per thread and do the per-hash work inside f: agents are not thread safe.
    template<typename function_t>
    void with_agent(function_t &&f) const {
        std::visit([&](auto const &ibf) {
            auto agent = ibf.membership_agent();
            f(agent);
        }, ibf_);
    }

    /*!\cond DEV
     * \brief Serialisation support. Version >= 4 files (preceded by format_magic) hold the layout and
     * a flag saying whether the IBF is compressed; older files hold a compressed IBF and no layout.
     */
    template<typename archive_t>
    void save(archive_t &archive) const {
        archive(window_size_);
        archive(kmer_size_);
        archive(max_fpr_);
        archive(summary_);
        archive(stats_);
        archive(layout_);
        const bool compressed = is_compressed();
        archive(compressed);
        std::visit([&](auto const &ibf) { archive(ibf); }, ibf_);
    }

    template<typename archive_t>
    void load(archive_t &archive, const bool has_magic) {
        try {
            archive(window_size_);
            archive(kmer_size_);
            archive(max_fpr_);
            archive(summary_);
            archive(stats_);
            bool compressed = true;
            if (has_magic) {
                archive(layout_);
                archive(compressed);
            } else {
                layout_ = Layout::identity(summary_.num_bins);
            }
            if (compressed) {
                ibf_compressed ibf;
                archive(ibf);
                ibf_ = std::move(ibf);
            } else {
                ibf_uncompressed ibf;
                archive(ibf);
                ibf_ = std::move(ibf);
            }
        }
            // GCOVR_EXCL_START
        catch (std::exception const &e) {
            PLOG_ERROR << "Cannot read index: " + std::string{e.what()};
            exit(1);
        }
    }
    //!\endcond
};

#endif // CHARON_INDEX_H
