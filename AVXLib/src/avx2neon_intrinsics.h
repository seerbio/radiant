#ifndef AVXLIB_AVX2NEON_INTRINSICS_H
#define AVXLIB_AVX2NEON_INTRINSICS_H

#include <arm_neon.h>
#include <cstdint>

/**
 * This file includes implementations of AVX2 intrinsics that are missing from AvxToNeon,
 * even though some are documented as being supported.
 */

// _mm256_loadu_ps: Load 8 floats (unaligned) from memory into __m256
static inline __m256 _mm256_loadu_ps(const float* mem_addr) {
    __m256 result;
    result.vect_f32[0] = vld1q_f32(mem_addr);       // Load first 4 floats
    result.vect_f32[1] = vld1q_f32(mem_addr + 4);   // Load next 4 floats
    return result;
}

// _mm256_storeu_ps: Store 8 floats (unaligned) from __m256 to memory
static inline void _mm256_storeu_ps(float* mem_addr, __m256 a) {
    vst1q_f32(mem_addr, a.vect_f32[0]);
    vst1q_f32(mem_addr + 4, a.vect_f32[1]);
}

// _mm256_unpacklo_epi16: Interleave lower 8 int16_t from each __m256i, zero-extend to int32_t
static inline __m256i _mm256_unpacklo_epi16(__m256i a, __m256i b) {
    __m256i result;
    result.vect_i128[0] = vreinterpretq_s64_s32(vmovl_s16(vget_low_s16(vreinterpretq_s16_s64(a.vect_i128[0]))));
    result.vect_i128[1] = vreinterpretq_s64_s32(vmovl_s16(vget_low_s16(vreinterpretq_s16_s64(b.vect_i128[0]))));
    return result;
}

// _mm256_unpackhi_epi16: Interleave upper 8 int16_t from each __m256i, zero-extend to int32_t
static inline __m256i _mm256_unpackhi_epi16(__m256i a, __m256i b) {
    __m256i result;
    result.vect_i128[0] = vreinterpretq_s64_s32(vmovl_s16(vget_high_s16(vreinterpretq_s16_s64(a.vect_i128[0]))));
    result.vect_i128[1] = vreinterpretq_s64_s32(vmovl_s16(vget_high_s16(vreinterpretq_s16_s64(b.vect_i128[0]))));
    return result;
}

#endif // AVXLIB_AVX2NEON_INTRINSICS_H
