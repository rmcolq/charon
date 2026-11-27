#ifndef CHARON_MD5_H
#define CHARON_MD5_H

#pragma once

// Small streaming MD5 (RFC 1321), used to check downloaded files against the checksums Zenodo publishes.
// Not for security purposes.

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>

class Md5 {
public:
    Md5() = default;

    void update(const void *data, std::size_t len) {
        const auto *bytes = static_cast<const uint8_t *>(data);
        std::size_t used = static_cast<std::size_t>(total_ % 64);
        total_ += len;
        if (used) {
            const std::size_t take = std::min(len, 64 - used);
            std::memcpy(buffer_.data() + used, bytes, take);
            bytes += take;
            len -= take;
            if (used + take < 64)
                return;
            transform(buffer_.data());
        }
        while (len >= 64) {
            transform(bytes);
            bytes += 64;
            len -= 64;
        }
        if (len)
            std::memcpy(buffer_.data(), bytes, len);
    }

    void update(std::string_view s) { update(s.data(), s.size()); }

    /// Finish and return the lower-case hex digest. The object must not be updated afterwards.
    std::string hexdigest() {
        const uint64_t bit_length = total_ * 8;
        static constexpr uint8_t pad_start = 0x80;
        static constexpr std::array<uint8_t, 64> zeros{};
        update(&pad_start, 1);
        const std::size_t used = static_cast<std::size_t>(total_ % 64);
        update(zeros.data(), used <= 56 ? 56 - used : 120 - used);
        std::array<uint8_t, 8> length_bytes{};
        for (int i = 0; i < 8; ++i)
            length_bytes[i] = static_cast<uint8_t>(bit_length >> (8 * i));
        update(length_bytes.data(), length_bytes.size());

        static constexpr char hex[] = "0123456789abcdef";
        std::string out;
        out.reserve(32);
        for (uint32_t word: state_)
            for (int i = 0; i < 4; ++i) {
                const auto byte = static_cast<uint8_t>(word >> (8 * i));
                out.push_back(hex[byte >> 4]);
                out.push_back(hex[byte & 0xf]);
            }
        return out;
    }

private:
    std::array<uint32_t, 4> state_{0x67452301u, 0xefcdab89u, 0x98badcfeu, 0x10325476u};
    std::array<uint8_t, 64> buffer_{};
    uint64_t total_{0};

    static constexpr uint32_t rotl(uint32_t x, uint32_t c) { return (x << c) | (x >> (32 - c)); }

    void transform(const uint8_t *block) {
        static constexpr std::array<uint32_t, 64> K{
                0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
                0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
                0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
                0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
                0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
                0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
                0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
                0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
        static constexpr std::array<uint32_t, 64> S{
                7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
                4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};

        std::array<uint32_t, 16> M{};
        for (int i = 0; i < 16; ++i)
            M[i] = static_cast<uint32_t>(block[4 * i]) | (static_cast<uint32_t>(block[4 * i + 1]) << 8) |
                   (static_cast<uint32_t>(block[4 * i + 2]) << 16) | (static_cast<uint32_t>(block[4 * i + 3]) << 24);

        uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        for (uint32_t i = 0; i < 64; ++i) {
            uint32_t f, g;
            if (i < 16) {
                f = (b & c) | (~b & d);
                g = i;
            } else if (i < 32) {
                f = (d & b) | (~d & c);
                g = (5 * i + 1) % 16;
            } else if (i < 48) {
                f = b ^ c ^ d;
                g = (3 * i + 5) % 16;
            } else {
                f = c ^ (b | ~d);
                g = (7 * i) % 16;
            }
            const uint32_t tmp = d;
            d = c;
            c = b;
            b = b + rotl(a + f + K[i] + M[g], S[i]);
            a = tmp;
        }
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
    }
};

#endif // CHARON_MD5_H
