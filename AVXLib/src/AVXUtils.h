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

	static Err convolveEightVecsWithKernelAVXFloat(
		const QVector<float> &kernel,
		size_t maxVecSize,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7
		);

	static Err findApexesEightVecs (
		size_t maxVecSize,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7,
		float *masterVectorApexes
		);

	static Err sumParallel (
		size_t maxVecSize,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7,
		float *sums
		);

	static Err findPeaksEightVecs (
		size_t maxVecSize,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7
		);

	static Err subtractArraysAVX2(
		float* array1,
		const float* array2,
		size_t size,
		bool zeroNegatives
		);

	static Err splitAVXUInt16to32(
		__m256i intShort,
		__m256i *int32Output1,
		__m256i *int32Output2
		);

	static bool isNByteAligned(
		const void* ptr,
		size_t byteSize
		) {
		return reinterpret_cast<uintptr_t>(ptr) % byteSize == 0;
	}

	static void replaceArrayValuesAVXGreaterThan(
		float threshold,
		float relacementValue,
		__m256 &arr
		);

	static void replaceArrayValuesAVXLessThan(
		float threshold,
		float relacementValue,
		__m256 &arr
		);

	static float cosineSimilarityAVX(
		__m256 vec1,
		__m256 vec2
		);

	static float dotProductAvx(
		const float* arrayA,
		const float* arrayB,
		size_t length
		);

	static float magnitudeAvx(
		const float* array,
		size_t length
		);

	static float cosineSimilarityAVX(
		const float* arrayA,
		const float* arrayB,
		size_t length
		);

	static void cosineSimilarityIntraAVXParallel(
		const float* arrRef,
		const float* array0,
		const float* array1,
		const float* array2,
		const float* array3,
		const float* array4,
		const float* array5,
		const float* array6,
		const float* array7,
		size_t length,
		float* cosineSimResultsAligned
		);

	static float maxFloat(__m256 vec);

	static bool isAllOnes(__m256 mask);

	static bool isAllZeros(__m256 vec);

	static void interleaveVectors(
		size_t size,
		size_t paddingSingleRegister,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7,
		float* resultVector
		);

	static void separateInterleavedVectors(
		const float* interleavedVectors,
		size_t size,
		size_t paddingSingleRegister,
		float* v0,
		float* v1,
		float* v2,
		float* v3,
		float* v4,
		float* v5,
		float* v6,
		float* v7
		);

	static __m256 log256(__m256 x);

	static float sum256(__m256 vec);

	static QVector<QPair<int, float>> findApexesAVX2(const float* array, size_t size);

};



#endif //AVXUTILS_H
