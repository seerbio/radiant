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

	static Err copyAVXFFloat(
		float* src,
		float* dst,
		size_t size

		) {

		ERR_INIT

		e = ErrorUtils::isByteAligned(dst, AVX2_ALIGNAS_SIZE);

		size_t i = 0;

		for (; i + AVX2_FLOAT_REGISTER_SIZE <= size; i += AVX2_FLOAT_REGISTER_SIZE) {
			__m256 v = _mm256_loadu_ps(src + i);
			_mm256_store_ps(dst + i, v);
		}

		if (i < size) {
			std::memcpy(dst + i, src + i, (size - i) * sizeof(float));
		}

		ERR_RETURN
	}



	static size_t calculateNextAlignedBlockSize(
		size_t arrSize,
		size_t byteSize
		);

};



#endif //AVXUTILS_H
