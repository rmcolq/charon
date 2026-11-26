#pragma once

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <plog/Log.h>
#include <cereal_unordered_dense.h>

struct InputSummary {
    uint8_t num_bins{0};
    std::vector<std::string> categories;
    std::vector<std::pair<std::string, uint8_t>> filepath_to_bin;
    ankerl::unordered_dense::map<uint8_t, std::string> bin_to_category;

    InputSummary() = default;
    InputSummary(InputSummary const&) = default;
    InputSummary(InputSummary&&) = default;
    InputSummary& operator=(InputSummary const&) = default;
    InputSummary& operator=(InputSummary&&) = default;
    ~InputSummary() = default;

    [[nodiscard]] uint8_t num_categories() const {
        return static_cast<uint8_t>(categories.size());
    }

    [[nodiscard]] uint8_t category_index(const std::string& category) const {
        for (size_t i = 0; i < categories.size(); ++i) {
            if (category == categories.at(i))
                return static_cast<uint8_t>(i);
        }
        return std::numeric_limits<uint8_t>::max();
    }

    [[nodiscard]] uint8_t host_category_index() const {
        const auto index1 = category_index("human");
        const auto index2 = category_index("host");
        const auto index = std::min(index1, index2);
        if (index == std::numeric_limits<uint8_t>::max()) {
            PLOG_ERROR << "Neither 'human' nor 'host' appear as categories in the index";
        }
        assert(index != std::numeric_limits<uint8_t>::max());
        return index;
    }

    // Fixed: was `index > categories.size()` which allowed index == size() through,
    // causing an out-of-bounds access in categories.at(index).
    [[nodiscard]] std::string category_name(const uint8_t index) const {
        if (index >= categories.size())
            return "";
        return categories.at(index);
    }

    // Replace each input file path with the shortest trailing part of it that no other input file shares:
    // the file name, or more of its directories when two files have the same name. A stored index then
    // doesn't record where its input files were, and can be moved to another machine.
    void shorten_filepaths() {
        std::vector<std::vector<std::string>> components;  // directories and file name, without the root
        components.reserve(filepath_to_bin.size());
        for (const auto& [path, bin] : filepath_to_bin) {
            std::vector<std::string> parts;
            for (const auto& part : std::filesystem::path{path}.relative_path())
                if (not part.empty() and part != ".")
                    parts.push_back(part.string());
            components.push_back(std::move(parts));
        }
        auto suffix = [&](const size_t i, const size_t length) {
            const auto& parts = components[i];
            std::filesystem::path result;
            for (size_t j = parts.size() - std::min(length, parts.size()); j < parts.size(); ++j)
                result /= parts[j];
            return result.generic_string();
        };

        std::vector<size_t> lengths(components.size(), 1);
        bool changed = true;
        while (changed) {
            changed = false;
            ankerl::unordered_dense::map<std::string, std::vector<size_t>> by_suffix;
            for (size_t i = 0; i < components.size(); ++i)
                by_suffix[suffix(i, lengths[i])].push_back(i);
            for (const auto& [name, indices] : by_suffix) {
                if (indices.size() < 2)
                    continue;
                for (const auto i : indices) {
                    if (lengths[i] < components[i].size()) {
                        ++lengths[i];
                        changed = true;
                    }
                }
            }
        }
        for (size_t i = 0; i < components.size(); ++i)
            filepath_to_bin[i].first = suffix(i, lengths[i]);
    }

    template<seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t& archive) {
        try {
            archive(num_bins);
            archive(categories);
            archive(filepath_to_bin);
            archive(bin_to_category);
        } catch (std::exception const& e) {
            PLOG_ERROR << "Cannot read input_summary: " + std::string{e.what()};
            exit(1);
        }
    }

    template<seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t& archive) {
        try {
            archive(num_bins);
            archive(categories);
            archive(filepath_to_bin);
            archive(bin_to_category);
        } catch (std::exception const& e) {
            PLOG_ERROR << "Cannot read input_summary: " + std::string{e.what()};
            exit(1);
        }
    }
};