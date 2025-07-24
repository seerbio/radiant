//
// Created by andrewnichols on 7/4/25.
//

#ifndef AVXUTILS_H
#define AVXUTILS_H

#include "AVXLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"

#include <cstdint>
#include <cstring>

#if defined(__aarch64__) || defined(__arm__) || defined(__ARM_NEON)
#warning "Building on ARM; will use avx2neon to translate intrinsics!"
#include <avx2neon.h>
#include "avx2neon_intrinsics.h" // Provide missing AVX2 intrinsics for ARM
#else
#include <immintrin.h>
#endif


using namespace Error;

class AVXLIB_EXPORTS AVXUtils {

public:

	static inline constexpr size_t AVX2_FLOAT_REGISTER_SIZE = 8;
	static inline constexpr size_t AVX2_INT16_REGISTER_SIZE = 16;
	static inline constexpr size_t AVX2_INT32_REGISTER_SIZE = 8;
	static inline constexpr size_t AVX2_ALIGNAS_SIZE = 32;

	static Err copyAVXFloatToAligned(
		float* src,
		float* dst,
		size_t size
		);

	static Err copyAVXIntToAligned(
		uint32_t *src,
		uint32_t *dst,
		size_t size
		);

	static size_t calculateNextAlignedBlockSize(
		size_t arrSize,
		size_t byteSize
		);

	static void printAVXFloat(const __m256 &avx);
	static void printAVXInt16(const __m256i &avx);
	static void printAVXInt32(const __m256i &avx);

	template<typename T>
	static void printMask(
		const T &mask,
		size_t maskSize
		) {

		QDebug debug = qDebug().nospace();;
		for (int j = 0; j < maskSize; j++) {
			debug << ((mask >> j) & 1) << ", ";
		}
	}

	static void printMask(const __m256& mask, size_t size);

	static Err convolveWithKernelAVXFloat(
		const QVector<float> &kernel,
		size_t size,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7
		);

	static Err splitAVXUInt16to32(
		__m256i intShort,
		__m256i *int32Output1,
		__m256i *int32Output2
		);


};



#endif //AVXUTILS_H
