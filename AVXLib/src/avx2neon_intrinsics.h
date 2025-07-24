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
    // Interleave lower 4 elements of each 128-bit lane
    int16x4_t a0 = vget_low_s16(a.vect_s16[0]);
    int16x4_t b0 = vget_low_s16(b.vect_s16[0]);
    int16x4_t a1 = vget_low_s16(a.vect_s16[1]);
    int16x4_t b1 = vget_low_s16(b.vect_s16[1]);
    // Interleave
    int16x4x2_t zipped0 = vzip_s16(a0, b0);
    int16x4x2_t zipped1 = vzip_s16(a1, b1);
    // Widen to int32x4_t
    result.vect_s32[0] = vmovl_s16(zipped0.val[0]);
    result.vect_s32[1] = vmovl_s16(zipped1.val[0]);
    return result;
}

// _mm256_unpackhi_epi16: Interleave upper 8 int16_t from each __m256i, zero-extend to int32_t
static inline __m256i _mm256_unpackhi_epi16(__m256i a, __m256i b) {
    __m256i result;
    // Interleave upper 4 elements of each 128-bit lane
    int16x4_t a0 = vget_high_s16(a.vect_s16[0]);
    int16x4_t b0 = vget_high_s16(b.vect_s16[0]);
    int16x4_t a1 = vget_high_s16(a.vect_s16[1]);
    int16x4_t b1 = vget_high_s16(b.vect_s16[1]);
    // Interleave
    int16x4x2_t zipped0 = vzip_s16(a0, b0);
    int16x4x2_t zipped1 = vzip_s16(a1, b1);
    // Widen to int32x4_t
    result.vect_s32[0] = vmovl_s16(zipped0.val[0]);
    result.vect_s32[1] = vmovl_s16(zipped1.val[0]);
    return result;
}

#endif // AVXLIB_AVX2NEON_INTRINSICS_H
