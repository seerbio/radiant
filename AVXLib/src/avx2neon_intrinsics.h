#ifndef AVXLIB_AVX2NEON_INTRINSICS_H
#define AVXLIB_AVX2NEON_INTRINSICS_H

#include <arm_neon.h>
#include <cstdint>

/**
 * This file includes implementations of AVX2 intrinsics that are missing from AvxToNeon,
 * even though some are documented as being supported.
 */

// _MM_SHUFFLE: Create 8-bit immediate for shuffle operations
#define _MM_SHUFFLE(z, y, x, w) (((z) << 6) | ((y) << 4) | ((x) << 2) | (w))

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

// _mm256_unpacklo_epi16: Interleave lower 8 uint16_t from each __m256i
static inline __m256i _mm256_unpacklo_epi16(__m256i a, __m256i b) {
    __m256i result;
    // Interleave lower 8 elements from each 128-bit lane
    uint16x8x2_t zipped0 = vzipq_u16(a.vect_u16[0], b.vect_u16[0]);
    uint16x8x2_t zipped1 = vzipq_u16(a.vect_u16[1], b.vect_u16[1]);
    result.vect_u16[0] = zipped0.val[0];
    result.vect_u16[1] = zipped1.val[0];
    return result;
}

// _mm256_unpackhi_epi16: Interleave upper 8 uint16_t from each __m256i
static inline __m256i _mm256_unpackhi_epi16(__m256i a, __m256i b) {
    __m256i result;
    // Interleave upper 8 elements from each 128-bit lane
    uint16x8x2_t zipped0 = vzipq_u16(a.vect_u16[0], b.vect_u16[0]);
    uint16x8x2_t zipped1 = vzipq_u16(a.vect_u16[1], b.vect_u16[1]);
    result.vect_u16[0] = zipped0.val[1];
    result.vect_u16[1] = zipped1.val[1];
    return result;
}

// _mm256_and_ps: Bitwise AND of packed single-precision (float) elements
static inline __m256 _mm256_and_ps(__m256 a, __m256 b) {
    __m256 result;
    result.vect_f32[0] = vreinterpretq_f32_u32(
        vandq_u32(vreinterpretq_u32_f32(a.vect_f32[0]), vreinterpretq_u32_f32(b.vect_f32[0])));
    result.vect_f32[1] = vreinterpretq_f32_u32(
        vandq_u32(vreinterpretq_u32_f32(a.vect_f32[1]), vreinterpretq_u32_f32(b.vect_f32[1])));
    return result;
}

// _mm256_permutevar8x32_ps: Shuffle 8 floats using variable indices
static inline __m256 _mm256_permutevar8x32_ps(__m256 a, __m256i idx) {
    __m256 result;
    // Extract float values into temporary array for indexing
    float temp_a[8];
    vst1q_f32(temp_a, a.vect_f32[0]);
    vst1q_f32(temp_a + 4, a.vect_f32[1]);

    // Extract index values into temporary array
    uint32_t temp_idx[8];
    vst1q_u32(temp_idx, idx.vect_u32[0]);
    vst1q_u32(temp_idx + 4, idx.vect_u32[1]);

    // Perform permutation with index masking (only lower 3 bits used)
    float temp_result[8];
    for (int i = 0; i < 8; i++) {
        temp_result[i] = temp_a[temp_idx[i] & 0x7];
    }

    // Load results back into NEON vectors
    result.vect_f32[0] = vld1q_f32(temp_result);
    result.vect_f32[1] = vld1q_f32(temp_result + 4);
    return result;
}

// _mm256_cvtss_f32: Extract the lowest single-precision float from __m256
static inline float _mm256_cvtss_f32(__m256 a) {
    return vgetq_lane_f32(a.vect_f32[0], 0);
}

// _mm256_max_ps: Element-wise maximum of packed single-precision (float) elements
static inline __m256 _mm256_max_ps(__m256 a, __m256 b) {
    __m256 result;
    result.vect_f32[0] = vmaxq_f32(a.vect_f32[0], b.vect_f32[0]);
    result.vect_f32[1] = vmaxq_f32(a.vect_f32[1], b.vect_f32[1]);
    return result;
}

// _mm256_castps_si256: Cast __m256 to __m256i (bitwise, no conversion)
static inline __m256i _mm256_castps_si256(__m256 a) {
    __m256i result;
    result.vect_u32[0] = vreinterpretq_u32_f32(a.vect_f32[0]);
    result.vect_u32[1] = vreinterpretq_u32_f32(a.vect_f32[1]);
    return result;
}

// _mm256_srli_epi32: Logical right shift of packed 32-bit integers by immediate
static inline __m256i _mm256_srli_epi32(__m256i a, int imm8) {
    __m256i result;
    result.vect_u32[0] = vshrq_n_u32(a.vect_u32[0], imm8);
    result.vect_u32[1] = vshrq_n_u32(a.vect_u32[1], imm8);
    return result;
}

// _mm256_fmadd_ps: Fused multiply-add of packed single-precision (float) elements (a * b + c)
static inline __m256 _mm256_fmadd_ps(__m256 a, __m256 b, __m256 c) {
    __m256 result;
    result.vect_f32[0] = vfmaq_f32(c.vect_f32[0], a.vect_f32[0], b.vect_f32[0]);
    result.vect_f32[1] = vfmaq_f32(c.vect_f32[1], a.vect_f32[1], b.vect_f32[1]);
    return result;
}

// _mm256_permute2f128_ps: Permute 128-bit lanes of packed single-precision floats
static inline __m256 _mm256_permute2f128_ps(__m256 a, __m256 b, int imm8) {
    __m256 result;

    // Extract control bits for each 128-bit lane
    int sel0 = imm8 & 0x3;        // bits [1:0] - select source for lower 128 bits
    int sel1 = (imm8 >> 4) & 0x3; // bits [5:4] - select source for upper 128 bits
    int zero0 = (imm8 >> 3) & 0x1; // bit 3 - zero lower 128 bits
    int zero1 = (imm8 >> 7) & 0x1; // bit 7 - zero upper 128 bits

    // Select lower 128-bit lane
    if (zero0) {
        result.vect_f32[0] = vdupq_n_f32(0.0f);
    } else {
        switch (sel0) {
            case 0: result.vect_f32[0] = a.vect_f32[0]; break; // a[127:0]
            case 1: result.vect_f32[0] = a.vect_f32[1]; break; // a[255:128]
            case 2: result.vect_f32[0] = b.vect_f32[0]; break; // b[127:0]
            case 3: result.vect_f32[0] = b.vect_f32[1]; break; // b[255:128]
        }
    }

    // Select upper 128-bit lane
    if (zero1) {
        result.vect_f32[1] = vdupq_n_f32(0.0f);
    } else {
        switch (sel1) {
            case 0: result.vect_f32[1] = a.vect_f32[0]; break; // a[127:0]
            case 1: result.vect_f32[1] = a.vect_f32[1]; break; // a[255:128]
            case 2: result.vect_f32[1] = b.vect_f32[0]; break; // b[127:0]
            case 3: result.vect_f32[1] = b.vect_f32[1]; break; // b[255:128]
        }
    }

    return result;
}

// _mm256_shuffle_ps: Shuffle single-precision floats within 128-bit lanes
static inline __m256 _mm256_shuffle_ps(__m256 a, __m256 b, int imm8) {
    __m256 result;

    // Extract shuffle control indices
    int w = imm8 & 0x3;        // bits [1:0] - select from b for element 0
    int x = (imm8 >> 2) & 0x3; // bits [3:2] - select from b for element 1
    int y = (imm8 >> 4) & 0x3; // bits [5:4] - select from a for element 2
    int z = (imm8 >> 6) & 0x3; // bits [7:6] - select from a for element 3

    // Shuffle lower 128-bit lane: {a[y], a[z], b[w], b[x]}
    float temp_low[4] = {
        vgetq_lane_f32(a.vect_f32[0], y),
        vgetq_lane_f32(a.vect_f32[0], z),
        vgetq_lane_f32(b.vect_f32[0], w),
        vgetq_lane_f32(b.vect_f32[0], x)
    };
    result.vect_f32[0] = vld1q_f32(temp_low);

    // Shuffle upper 128-bit lane: {a[y], a[z], b[w], b[x]}
    float temp_high[4] = {
        vgetq_lane_f32(a.vect_f32[1], y),
        vgetq_lane_f32(a.vect_f32[1], z),
        vgetq_lane_f32(b.vect_f32[1], w),
        vgetq_lane_f32(b.vect_f32[1], x)
    };
    result.vect_f32[1] = vld1q_f32(temp_high);

    return result;
}

// _mm_cvtss_f32: Extract the lowest single-precision float from __m128
static inline float _mm_cvtss_f32(__m128 a) {
    return vgetq_lane_f32(a, 0);
}

// _mm256_sqrt_ps: Square root of packed single-precision (float) elements
static inline __m256 _mm256_sqrt_ps(__m256 a) {
    __m256 result;
    result.vect_f32[0] = vsqrtq_f32(a.vect_f32[0]);
    result.vect_f32[1] = vsqrtq_f32(a.vect_f32[1]);
    return result;
}

// _mm_add_ps: Add packed single-precision (float) elements
static inline __m128 _mm_add_ps(__m128 a, __m128 b) {
    return vaddq_f32(a, b);
}

// _mm_hadd_ps: Horizontal add of packed single-precision (float) elements
static inline __m128 _mm_hadd_ps(__m128 a, __m128 b) {
    float32x4_t a_temp = vpaddq_f32(a, a);    // {a0+a1, a2+a3, a0+a1, a2+a3}
    float32x4_t b_temp = vpaddq_f32(b, b);    // {b0+b1, b2+b3, b0+b1, b2+b3}
    return vcombine_f32(vget_low_f32(a_temp), vget_low_f32(b_temp));  // {a0+a1, a2+a3, b0+b1, b2+b3}
}

// _mm_max_ps: Element-wise maximum of packed single-precision (float) elements
static inline __m128 _mm_max_ps(__m128 a, __m128 b) {
    return vmaxq_f32(a, b);
}

// _mm_movehl_ps: Move high two single-precision floats from a to low two positions of result
static inline __m128 _mm_movehl_ps(__m128 a, __m128 b) {
    return vcombine_f32(vget_high_f32(b), vget_high_f32(a));
}

// _mm_shuffle_ps: Shuffle single-precision floats using immediate control
static inline __m128 _mm_shuffle_ps(__m128 a, __m128 b, int imm8) {
    // Extract shuffle control indices
    int w = imm8 & 0x3;        // bits [1:0] - select from a for element 0
    int x = (imm8 >> 2) & 0x3; // bits [3:2] - select from a for element 1
    int y = (imm8 >> 4) & 0x3; // bits [5:4] - select from b for element 2
    int z = (imm8 >> 6) & 0x3; // bits [7:6] - select from b for element 3

    // Result: {a[w], a[x], b[y], b[z]}
    float temp[4] = {
        vgetq_lane_f32(a, w),
        vgetq_lane_f32(a, x),
        vgetq_lane_f32(b, y),
        vgetq_lane_f32(b, z)
    };
    return vld1q_f32(temp);
}

// _mm256_testc_si256: Test if all bits in (a AND b) are clear (return 1 if all zero, 0 otherwise)
static inline int _mm256_testc_si256(__m256i a, __m256i b) {
    // Perform bitwise AND on both 128-bit lanes
    uint32x4_t and_low = vandq_u32(a.vect_u32[0], b.vect_u32[0]);
    uint32x4_t and_high = vandq_u32(a.vect_u32[1], b.vect_u32[1]);

    // Check if all bits are zero in both lanes
    uint64x2_t or_low = vreinterpretq_u64_u32(and_low);
    uint64x2_t or_high = vreinterpretq_u64_u32(and_high);

    // Combine lanes and test for all zeros
    uint64x2_t combined = vorrq_u64(or_low, or_high);
    uint64_t result = vgetq_lane_u64(combined, 0) | vgetq_lane_u64(combined, 1);

    return (result == 0) ? 1 : 0;
}

#endif // AVXLIB_AVX2NEON_INTRINSICS_H
