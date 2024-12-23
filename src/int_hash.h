/* C++ implementation of the Robert Jenkins' 64-bit Mix Function with Seed */

#ifndef __INTHASH_H__
#define __INTHASH_H__

#include <cstdint>

//A series of arithmetic (+, -) and bitwise operations (<<, >>, ^) are applied to mix the bits thoroughly.
inline uint64_t jenkins_hash_seed(uint64_t x, uint64_t seed) {
    x ^= seed;
    x = (x + 0x7ed55d16ULL) + (x << 12);
    x = (x ^ 0xc761c23cULL) ^ (x >> 19);
    x = (x + 0x165667b1ULL) + (x << 5);
    x = (x + 0xd3a2646cULL) ^ (x << 9);
    x = (x + 0xfd7046c5ULL) + (x << 3);
    x = (x ^ 0xb55a4f09ULL) ^ (x >> 16);
    return x;
}

// The unhash function reverses the operations applied in the hash function in the exact reverse order.
inline uint64_t jenkins_unhash_seed(uint64_t h, uint64_t seed) {
    h = (h ^ 0xb55a4f09ULL) ^ (h >> 16);
    h = (h - 0xfd7046c5ULL) - (h << 3);
    h = (h ^ 0xd3a2646cULL) ^ (h >> 9);
    h = (h - 0x165667b1ULL) - (h << 5);
    h = (h ^ 0xc761c23cULL) ^ (h >> 19);
    h = (h - 0x7ed55d16ULL) - (h << 12);
    return h ^ seed;
}

inline uint64_t
murmur3_fmix(uint64_t x, uint64_t seed) noexcept {
	x ^= seed; 
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccd;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53;
    x ^= x >> 33;
    return x;
}

inline uint32_t
murmur3_fmix(uint32_t x, uint32_t seed) noexcept {
	x ^= seed;
    x ^= x >> 16;
    x *= 0x85ebca6b;
    x ^= x >> 13;
    x *= 0xc2b2ae35;
    x ^= x >> 16;
    return x;
}
#endif
