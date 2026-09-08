#ifndef CHARON_CEREAL_UNORDERED_DENSE_H
#define CHARON_CEREAL_UNORDERED_DENSE_H

#pragma once

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <vector>
#include <utility>

// This header must be included before this one in files that use ankerl maps
// We don't include it here to avoid circular dependencies

namespace cereal {

// Serialization support for ankerl::unordered_dense::map
// Specialization must be in cereal namespace
template <class Archive, class Key, class T, class Hash, class KeyEqual, class Allocator, class BucketType>
void CEREAL_SAVE_FUNCTION_NAME(Archive& ar, ankerl::unordered_dense::map<Key, T, Hash, KeyEqual, Allocator, BucketType> const& map) {
    // Save as vector of pairs
    std::vector<std::pair<Key, T>> items(map.begin(), map.end());
    ar(items);
}

template <class Archive, class Key, class T, class Hash, class KeyEqual, class Allocator, class BucketType>
void CEREAL_LOAD_FUNCTION_NAME(Archive& ar, ankerl::unordered_dense::map<Key, T, Hash, KeyEqual, Allocator, BucketType>& map) {
    // Load from vector of pairs
    std::vector<std::pair<Key, T>> items;
    ar(items);
    map.clear();
    map.reserve(items.size());
    for (auto& item : items) {
        map.emplace(std::move(item.first), std::move(item.second));
    }
}

// Serialization support for ankerl::unordered_dense::set
template <class Archive, class Key, class Hash, class KeyEqual, class Allocator, class BucketType>
void CEREAL_SAVE_FUNCTION_NAME(Archive& ar, ankerl::unordered_dense::set<Key, Hash, KeyEqual, Allocator, BucketType> const& set) {
    // Save as vector
    std::vector<Key> items(set.begin(), set.end());
    ar(items);
}

template <class Archive, class Key, class Hash, class KeyEqual, class Allocator, class BucketType>
void CEREAL_LOAD_FUNCTION_NAME(Archive& ar, ankerl::unordered_dense::set<Key, Hash, KeyEqual, Allocator, BucketType>& set) {
    // Load from vector
    std::vector<Key> items;
    ar(items);
    set.clear();
    set.reserve(items.size());
    for (auto& item : items) {
        set.emplace(std::move(item));
    }
}

} // namespace cereal

#endif // CHARON_CEREAL_UNORDERED_DENSE_H
