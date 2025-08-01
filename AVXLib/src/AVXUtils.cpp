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

void AVXUtils::interleaveVectors(
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
	) {

	for (int i = 0; i < size; i++) {
		__m256 parallelVec = _mm256_set_ps(v7[i], v6[i], v5[i], v4[i], v3[i], v2[i], v1[i], v0[i]);
		const size_t insertIndex = (paddingSingleRegister * AVX2_FLOAT_REGISTER_SIZE) + (i * AVX2_FLOAT_REGISTER_SIZE);
		_mm256_storeu_ps(resultVector + insertIndex, parallelVec);
	}
}

void AVXUtils::separateInterleavedVectors(
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
	) {

	const size_t registerCount = size / AVX2_FLOAT_REGISTER_SIZE;

	for (int i = static_cast<int>(paddingSingleRegister); i < registerCount  - paddingSingleRegister; i++) {
		v0[i - paddingSingleRegister] = interleavedVectors[i * AVX2_FLOAT_REGISTER_SIZE];
		v1[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+1];
		v2[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+2];
		v3[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+3];
		v4[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+4];
		v5[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+5];
		v6[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+6];
		v7[i - paddingSingleRegister] = interleavedVectors[(i * AVX2_FLOAT_REGISTER_SIZE)+7];
	}
}

Err AVXUtils::convolveEightVecsWithKernelAVXFloat(
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
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(kernel); ree;

	const int kernelSize = kernel.size();

	const size_t paddingSingleSide = (kernelSize - 1) / 2;
	const size_t paddingRegister = paddingSingleSide * AVX2_FLOAT_REGISTER_SIZE;
	const size_t masterVectorSize = (maxVecSize * AVX2_FLOAT_REGISTER_SIZE) + (paddingRegister * 2);
	const size_t resultsSize = maxVecSize * AVX2_FLOAT_REGISTER_SIZE;

	const size_t masterVectorSizeAlignas = calculateNextAlignedBlockSize(
		masterVectorSize,
		AVX2_ALIGNAS_SIZE
		);

	alignas(AVX2_ALIGNAS_SIZE) float masterVector[masterVectorSizeAlignas] = {0};
	alignas(AVX2_ALIGNAS_SIZE) float result[masterVectorSizeAlignas] = {0};

	interleaveVectors(
		 maxVecSize,
		 paddingSingleSide,
		 v0,
		 v1,
		 v2,
		 v3,
		 v4,
		 v5,
		 v6,
		 v7,
		 masterVector
		);

	__m256 previousRegisters[kernelSize];
	for (int i = 0; i < kernelSize; i++) {
		previousRegisters[i] = _mm256_setzero_ps();
	}

	__m256 filterRegisters[kernelSize];
	for (int i = 0; i < kernelSize; i++) {
		filterRegisters[i] = _mm256_set1_ps(kernel[i]);
	}

	for (int i = paddingSingleSide; i < maxVecSize  + paddingSingleSide; i++) {

		const size_t readIndex = i * AVX2_FLOAT_REGISTER_SIZE;
		const __m256 registerCurrent = _mm256_load_ps(masterVector + readIndex);

		for (int j = 0; j < kernelSize - 1; j++) {
			previousRegisters[j] = previousRegisters[j+1];
		}
		previousRegisters[kernelSize - 1] = registerCurrent;

		__m256 vSum = _mm256_setzero_ps();
		for (size_t j = 0; j < kernelSize; j++) {
			const __m256 multVal = _mm256_mul_ps(previousRegisters[j], filterRegisters[j]);
			vSum = _mm256_add_ps(vSum, multVal);
		}

		const size_t writeIndex = (i - (paddingSingleSide * 2)) * AVX2_FLOAT_REGISTER_SIZE;
		_mm256_store_ps(result + writeIndex, vSum);
	}

	separateInterleavedVectors(
		result,
		resultsSize,
		0,
		v0,
		v1,
		v2,
		v3,
		v4,
		v5,
		v6,
		v7
		);

	ERR_RETURN
}

Err AVXUtils::findApexesEightVecs(
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
	) {

	ERR_INIT

	const size_t masterVectorSize = maxVecSize * AVX2_FLOAT_REGISTER_SIZE;

	const size_t masterVectorSizeAlignas = calculateNextAlignedBlockSize(
		masterVectorSize,
		AVX2_ALIGNAS_SIZE
		);

	alignas(AVX2_ALIGNAS_SIZE) float masterVector[masterVectorSizeAlignas] = {0};

	interleaveVectors(
		 maxVecSize,
		 0,
		 v0,
		 v1,
		 v2,
		 v3,
		 v4,
		 v5,
		 v6,
		 v7,
		 masterVector
		);

	__m256 zeroVec = _mm256_setzero_ps();

	__m256 previous = _mm256_load_ps(masterVector + 0);
	__m256 current = _mm256_load_ps(masterVector + AVX2_FLOAT_REGISTER_SIZE);
	for (int i = 2; i < maxVecSize; i++) {
		__m256 next = _mm256_load_ps(masterVector + (i * AVX2_FLOAT_REGISTER_SIZE));

		__m256 gtPrevious = _mm256_cmp_ps(current, previous, _CMP_GT_OQ);
		__m256 gtNext = _mm256_cmp_ps(current, next, _CMP_GT_OQ);
		__m256 apexMask = _mm256_and_ps(gtPrevious, gtNext);

		__m256 apexValues = _mm256_blendv_ps(zeroVec, current, apexMask);

		_mm256_store_ps(masterVectorApexes + ((i -1) * AVX2_FLOAT_REGISTER_SIZE), apexValues);

		previous = current;
		current = next;
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

void AVXUtils::replaceArrayValuesAVXGreaterThan(float threshold, float relacementValue, __m256 &arr) {
	const __m256 thresholds = _mm256_set1_ps(threshold);
	const __m256 mask = _mm256_cmp_ps(arr, thresholds, _CMP_GT_OQ);
	const __m256 replaceVals = _mm256_set1_ps(relacementValue);
	arr = _mm256_and_ps(mask, replaceVals);
}

void AVXUtils::replaceArrayValuesAVXLessThan(float threshold, float relacementValue, __m256 &arr) {
	const __m256 thresholds = _mm256_set1_ps(threshold);
	const __m256 replaceVals = _mm256_set1_ps(relacementValue);
	__m256 mask = _mm256_cmp_ps(arr, thresholds, _CMP_LT_OQ);
	arr = _mm256_blendv_ps(arr, replaceVals, mask);
}

namespace {
	float dotProduct256(__m256 a, __m256 b) {
		__m256 mul = _mm256_mul_ps(a, b);
		__m128 low = _mm256_castps256_ps128(mul);
		__m128 high = _mm256_extractf128_ps(mul, 1);
		__m128 sum128 = _mm_add_ps(low, high);
		sum128 = _mm_hadd_ps(sum128, sum128);
		sum128 = _mm_hadd_ps(sum128, sum128);
		return _mm_cvtss_f32(sum128);
	}
}
float AVXUtils::cosineSimilarityAVX(
	__m256 vec1,
	__m256 vec2
	) {

	const float dotProduct = dotProduct256(vec1, vec2);
	const float mag1 = std::sqrt(dotProduct256(vec1, vec1));
	const float mag2 = std::sqrt(dotProduct256(vec2, vec2));

	const float cosineSimilarity = dotProduct / (mag1 * mag2);
	return std::isnan(cosineSimilarity) ? 0.0f : cosineSimilarity;
}

float AVXUtils::dotProductAvx(const float *arrayA, const float *arrayB, size_t length) {

	__m256 sumVec = _mm256_setzero_ps();

	for (size_t i = 0; i < length; i += 8) {
		__m256 vecA = _mm256_load_ps(&arrayA[i]);
		__m256 vecB = _mm256_load_ps(&arrayB[i]);
		__m256 mul = _mm256_mul_ps(vecA, vecB);
		sumVec = _mm256_add_ps(sumVec, mul);
	}

	__m128 low = _mm256_castps256_ps128(sumVec);
	__m128 high = _mm256_extractf128_ps(sumVec, 1);
	__m128 sum128 = _mm_add_ps(low, high);
	sum128 = _mm_hadd_ps(sum128, sum128);
	sum128 = _mm_hadd_ps(sum128, sum128);

	return _mm_cvtss_f32(sum128);
}

float AVXUtils::magnitudeAvx(
	const float *array,
	size_t length
	) {
	return std::sqrt(dotProductAvx(array, array, length));
}

float AVXUtils::cosineSimilarityAVX(
	const float *arrayA,
	const float *arrayB,
	size_t length
	) {

	float dot = dotProductAvx(arrayA, arrayB, length);
	float magA = magnitudeAvx(arrayA, length);
	float magB = magnitudeAvx(arrayB, length);

	const float cosineSimilarity = dot / (magA * magB);
	return std::isnan(cosineSimilarity) ? 0.0f : cosineSimilarity;
}

float AVXUtils::maxFloat(__m256 vec)  {
	__m128 high = _mm256_extractf128_ps(vec, 1);
	__m128 low = _mm256_castps256_ps128(vec);

	__m128 max128 = _mm_max_ps(low, high);

	max128 = _mm_max_ps(max128, _mm_movehl_ps(max128, max128));
	max128 = _mm_max_ps(max128, _mm_shuffle_ps(max128, max128, 0x1));

	return _mm_cvtss_f32(max128);
}

bool AVXUtils::isAllOnes(__m256 mask) {
	__m256i intMask = _mm256_castps_si256(mask);
	__m256i allOnes = _mm256_set1_epi32(-1);
	return _mm256_testc_si256(intMask, allOnes);
}

bool AVXUtils::isAllZeros(__m256 vec) {
	__m256 zero_vec = _mm256_setzero_ps();
	__m256 cmp = _mm256_cmp_ps(vec, zero_vec, _CMP_EQ_OQ); // _CMP_EQ_OQ: Equal, ordered, quiet NaN handling
	int mask = _mm256_movemask_ps(cmp);
	return (mask == 0xFF);
}








