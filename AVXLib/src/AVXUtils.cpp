//
// Created by andrewnichols on 7/4/25.
//

#include "AVXUtils.h"

#include "GlobalSettings.h"


Err AVXUtils::copyAVXFloatToAligned(float *src, float *dst, size_t size) {

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

Err AVXUtils::copyAVXIntToAligned(uint32_t *src, uint32_t *dst, size_t size) {

	ERR_INIT

	e = ErrorUtils::isByteAligned(dst, AVX2_ALIGNAS_SIZE);

	size_t i = 0;

	for (; i + AVX2_INT32_REGISTER_SIZE <= size; i += AVX2_INT32_REGISTER_SIZE) {
		__m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + i));
		_mm256_store_si256(reinterpret_cast<__m256i*>(dst + i), v);
	}

	if (i < size) {
		std::memcpy(dst + i, src + i, (size - i) * sizeof(float));
	}

	ERR_RETURN
}

size_t AVXUtils::calculateNextAlignedBlockSize(
	size_t arrSize,
	size_t byteSize
	) {

	while (arrSize % byteSize != 0) {
		arrSize++;
	}

	return arrSize;
}

void AVXUtils::printAVXFloat(const __m256 &avx) {

	alignas(AVX2_ALIGNAS_SIZE) float values[AVX2_FLOAT_REGISTER_SIZE];
	_mm256_storeu_ps(values, avx);

	QDebug debug = qDebug().nospace();
	for (int j = 0; j < AVX2_FLOAT_REGISTER_SIZE; j++) {
		debug << values[j] << ", ";
	}
}

void AVXUtils::printAVXInt16(const __m256i &avx) {

	alignas(AVX2_ALIGNAS_SIZE) int16_t values[AVX2_INT16_REGISTER_SIZE];
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(values), avx);

	QDebug debug = qDebug().nospace();
	for (int j = 0; j < AVX2_INT16_REGISTER_SIZE; j++) {
		debug << values[j] << ", ";
	}
}

void AVXUtils::printAVXInt32(const __m256i &avx) {

	alignas(AVX2_ALIGNAS_SIZE) int32_t values[AVX2_INT32_REGISTER_SIZE];
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(values), avx);

	QDebug debug = qDebug().nospace();
	for (int j = 0; j < AVX2_INT32_REGISTER_SIZE; j++) {
		debug << values[j] << ", ";
	}
}

void AVXUtils::printMask(const __m256 &mask, size_t size) {
	int intMask = _mm256_movemask_ps(mask);

	QDebug debug = qDebug().nospace();;
	for (size_t j = 0; j < size; ++j) {
		debug << ((intMask >> j) & 1) << ", ";
	}
}

namespace {


}//namespace
Err AVXUtils::convolveWithKernelAVXFloat(
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
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(kernel); ree;

	const size_t padding = kernel.size() - 1;
	const size_t paddingSingle = padding / 2;
	const size_t paddingSingleRegister = paddingSingle * AVX2_FLOAT_REGISTER_SIZE;
	const size_t masterVectorSize = (size * AVX2_FLOAT_REGISTER_SIZE) + padding;

	const size_t masterVectorSizeAlignas = calculateNextAlignedBlockSize(
		masterVectorSize,
		AVX2_ALIGNAS_SIZE
		);

	alignas(AVX2_ALIGNAS_SIZE) float masterVector[masterVectorSizeAlignas] = {0};
	alignas(AVX2_ALIGNAS_SIZE) float result[masterVectorSize] = {0};

	for (int i = 0; i < size; i++) {
		__m256 parallelVec = _mm256_set_ps(v7[i], v6[i], v5[i], v4[i], v3[i], v2[i], v1[i], v0[i]);
		const size_t insertIndex = paddingSingleRegister + (i * AVX2_FLOAT_REGISTER_SIZE);
		_mm256_store_ps(masterVector + insertIndex, parallelVec);
	}

	for (int i = 0; i < (size * AVX2_FLOAT_REGISTER_SIZE); i += AVX2_FLOAT_REGISTER_SIZE) {

		__m256 vSum = _mm256_setzero_ps();

		for (int j = 0; j < kernel.size(); j++) {
			const size_t readIndex = i + (j * AVX2_FLOAT_REGISTER_SIZE);
			const __m256 reg = _mm256_load_ps(masterVector + readIndex);
			const __m256 parallelVec = _mm256_set1_ps(kernel[j]);
			const __m256 resultLocal = _mm256_mul_ps(reg, parallelVec);
			vSum = _mm256_add_ps(resultLocal, vSum);
		}
		_mm256_store_ps(result + i, vSum);
	}

	for (int i = 0; i < size; i++) {
		v0[i] = result[i];
		v1[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+1];
		v2[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+2];
		v3[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+3];
		v4[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+4];
		v5[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+5];
		v6[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+6];
		v7[i] = result[(i * AVX2_FLOAT_REGISTER_SIZE)+7];
	}

	ERR_RETURN
}

Err AVXUtils::subtractArraysAVX2(
	float *source,
	const float *subtractor,
	size_t size,
	bool zeroNegatives
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(isNByteAligned(source, AVX2_ALIGNAS_SIZE)); ree;
	e = ErrorUtils::isTrue(isNByteAligned(subtractor, AVX2_ALIGNAS_SIZE)); ree;

	const size_t avx2Iterations = size / AVX2_FLOAT_REGISTER_SIZE;

	for (size_t i = 0; i < avx2Iterations * AVX2_FLOAT_REGISTER_SIZE; i += AVX2_FLOAT_REGISTER_SIZE) {
		const __m256 a = _mm256_load_ps(source + i);
		const __m256 b = _mm256_load_ps(subtractor + i);
		__m256 resultVector = _mm256_sub_ps(a, b);

		if (zeroNegatives) {
			constexpr float thresholdValue = 0.0;
			constexpr float replaceValue = 0.0;
			replaceArrayValuesAVXLessThan(thresholdValue, replaceValue, resultVector);
		}

		_mm256_store_ps(source + i, resultVector);
	}

	for (size_t i = avx2Iterations * AVX2_FLOAT_REGISTER_SIZE; i < size; ++i) {
		source[i] = zeroNegatives
				  ? std::max(source[i] - subtractor[i], 0.0f)
				  : source[i] + subtractor[i];
	}

	ERR_RETURN
}

Err AVXUtils::splitAVXUInt16to32(
	__m256i intShort,
	__m256i *int32Output1,
	__m256i *int32Output2
	) {

	ERR_INIT

	*int32Output1 = _mm256_unpacklo_epi16(intShort, _mm256_setzero_si256());
	*int32Output2 = _mm256_unpackhi_epi16(intShort, _mm256_setzero_si256());

	ERR_RETURN
}








