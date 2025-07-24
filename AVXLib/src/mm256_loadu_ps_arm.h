#ifndef AVXLIB_MM256_LOADU_PS_H
#define AVXLIB_MM256_LOADU_PS_H

#include <arm_neon.h>
#include <cstdint>

// Emulate _mm256_loadu_ps for ARM NEON
// Loads 8 floats (unaligned) from memory into a float32x4x2_t, then combines into __m256 (float[8])
// This matches the approach used in avx2neon for other unaligned loads

static inline float32x4x2_t _mm256_loadu_ps(const float* mem_addr) {
    float32x4x2_t result;
    result.val[0] = vld1q_f32(mem_addr);       // Load first 4 floats
    result.val[1] = vld1q_f32(mem_addr + 4);   // Load next 4 floats
    return result;
}

// If you need to use as __m256, you may need a union or reinterpret_cast depending on your AVX2NEON setup

#endif // AVXLIB_MM256_LOADU_PS_H

