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
#include <immintrin.h>


using namespace Error;

class AVXLIB_EXPORTS AVXUtils {

public:

	static inline constexpr size_t AVX2_FLOAT_REGISTER_SIZE = 8;
	static inline constexpr size_t AVX2_ALIGNAS_SIZE = 32;

	static Err copyAVXFloatToAligned(
		float* src,
		float* dst,
		size_t size
		);

	static size_t calculateNextAlignedBlockSize(
		size_t arrSize,
		size_t byteSize
		);

	static void printAVXFloat(const __m256 &avx);

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
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7,
		size_t size
		);

};



#endif //AVXUTILS_H
