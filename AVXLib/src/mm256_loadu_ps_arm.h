#ifndef AVXLIB_MM256_LOADU_PS_H
#define AVXLIB_MM256_LOADU_PS_H

#include <arm_neon.h>
#include <cstdint>
#include "avxintrin.h" // for __m256 definition

// Emulate _mm256_loadu_ps for ARM NEON, matching AvxToNeon __m256 layout
static inline __m256 _mm256_loadu_ps(const float* mem_addr) {
    __m256 result;
    result.vect_f32[0] = vld1q_f32(mem_addr);       // Load first 4 floats
    result.vect_f32[1] = vld1q_f32(mem_addr + 4);   // Load next 4 floats
    return result;
}

#endif // AVXLIB_MM256_LOADU_PS_H
