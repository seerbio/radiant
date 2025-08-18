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
	constexpr int avx2 = AVX2_FLOAT_REGISTER_SIZE * 2;
	constexpr int avx3 = AVX2_FLOAT_REGISTER_SIZE * 3;

	for (size_t i = 0; i + 3 < size; i += 4) {
		__m256 p0 = _mm256_set_ps(v7[i],   v6[i],   v5[i],   v4[i],   v3[i],   v2[i],   v1[i],   v0[i]);
		__m256 p1 = _mm256_set_ps(v7[i+1], v6[i+1], v5[i+1], v4[i+1], v3[i+1], v2[i+1], v1[i+1], v0[i+1]);
		__m256 p2 = _mm256_set_ps(v7[i+2], v6[i+2], v5[i+2], v4[i+2], v3[i+2], v2[i+2], v1[i+2], v0[i+2]);
		__m256 p3 = _mm256_set_ps(v7[i+3], v6[i+3], v5[i+3], v4[i+3], v3[i+3], v2[i+3], v1[i+3], v0[i+3]);

		size_t base = (paddingSingleRegister + i) * AVX2_FLOAT_REGISTER_SIZE;
		_mm256_store_ps(resultVector + base,  p0);
		_mm256_store_ps(resultVector + base + AVX2_FLOAT_REGISTER_SIZE,  p1);
		_mm256_store_ps(resultVector + base + avx2, p2);
		_mm256_store_ps(resultVector + base + avx3, p3);
	}
	for (size_t i = (size & ~3); i < size; i++) {
		__m256 parallelVec = _mm256_set_ps(v7[i], v6[i], v5[i], v4[i], v3[i], v2[i], v1[i], v0[i]);
		const size_t insertIndex = (paddingSingleRegister + i) * AVX2_FLOAT_REGISTER_SIZE;
		_mm256_store_ps(resultVector + insertIndex, parallelVec);
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

	for (size_t i = paddingSingleRegister; i < registerCount - paddingSingleRegister; ++i) {
		__m256 vec = _mm256_load_ps(&interleavedVectors[i * AVX2_FLOAT_REGISTER_SIZE]);

		float* dst[] = { v0, v1, v2, v3, v4, v5, v6, v7 };
		for (int lane = 0; lane < AVX2_FLOAT_REGISTER_SIZE; ++lane) {
			dst[lane][i - paddingSingleRegister] = _mm256_cvtss_f32(_mm256_permutevar8x32_ps(vec, _mm256_set1_epi32(lane)));
		}
	}
}

__m256 AVXUtils::log256(__m256 x) {
	__m256 one = _mm256_set1_ps(1.0f);

	x = _mm256_max_ps(x, _mm256_set1_ps(1.17549e-38f)); // Avoid log(0)

	__m256i ix = _mm256_castps_si256(x);
	__m256i emm0 = _mm256_srli_epi32(ix, 23);

	// Keep only the fractional part
	ix = _mm256_and_si256(ix, _mm256_set1_epi32(0x7FFFFF));
	ix = _mm256_or_si256(ix, _mm256_set1_epi32(0x3f000000));
	__m256 f = _mm256_castsi256_ps(ix);

	// Exponent
	emm0 = _mm256_sub_epi32(emm0, _mm256_set1_epi32(127));
	__m256 e = _mm256_cvtepi32_ps(emm0);

	f = _mm256_sub_ps(f, one);

	// Approximate polynomial
	__m256 p = _mm256_set1_ps(7.0376836292e-2f);
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(-1.1514610310e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(1.1676998740e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(-1.2420140846e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(+1.4249322787e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(-1.6668057665e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(+2.0000714765e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(-2.4999993993e-1f));
	p = _mm256_fmadd_ps(p, f, _mm256_set1_ps(+3.3333331174e-1f));

	__m256 log1pf = _mm256_mul_ps(p, f);

	return _mm256_add_ps(_mm256_mul_ps(e, _mm256_set1_ps(0.69314718056f)), log1pf);
}

float AVXUtils::sum256(__m256 vec)  {
	__m256 shuff1 = _mm256_permute2f128_ps(vec, vec, 1);
	__m256 sum1 = _mm256_add_ps(vec, shuff1);

	__m256 shuff2 = _mm256_shuffle_ps(sum1, sum1, _MM_SHUFFLE(2, 3, 0, 1));
	__m256 sum2 = _mm256_add_ps(sum1, shuff2);

	__m256 shuff3 = _mm256_shuffle_ps(sum2, sum2, _MM_SHUFFLE(1, 0, 3, 2));
	__m256 sum3 = _mm256_add_ps(sum2, shuff3);

	return _mm_cvtss_f32(_mm256_castps256_ps128(sum3));
}

QVector<QPair<int, float>> AVXUtils::findApexesAVX2(const float *array, size_t size)  {

	QVector<QPair<int, float>> maxima;

	// if (size < 3) {
	// 	return maxima; // No apexes in array with fewer than 3 elements
	// }
	//
	// for (int i = 1; i < size - 1; ++i) {
	// 	if (array[i] > array[i - 1] && array[i] > array[i + 1]) {
	// 		maxima.push_back({i, array[i]});
	// 	}
	// }
	//
	// return maxima;

	if (size < 3) {
		return maxima;
	}

	int i = 1;
	for (; i + AVX2_FLOAT_REGISTER_SIZE < size - 1; i += 8) {

		const __m256 left = _mm256_loadu_ps(array + i - 1);
		const __m256 center = _mm256_loadu_ps(array + i);
		const __m256 right = _mm256_loadu_ps(array + i + 1);

		const __m256 greaterThanLeft = _mm256_cmp_ps(center, left, _CMP_GT_OS);
		const __m256 greaterThanRight = _mm256_cmp_ps(center, right, _CMP_GT_OS);
		const __m256 maximaMask = _mm256_and_ps(greaterThanLeft, greaterThanRight);

		const int mask = _mm256_movemask_ps(maximaMask);

		for (int j = 0; j < 8; ++j) {
			if (mask & (1 << j)) {
				maxima.push_back({i + j, array[i + j]});
			}
		}
	}

	for (; i < size - 1; ++i) {
		if (array[i] > array[i - 1] && array[i] > array[i + 1]) { // Local maxima
			maxima.push_back({i, array[i]});
		}
	}

    return maxima;
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

Err AVXUtils::findPeaksEightVecs(
	size_t maxVecSize,
	float *v0,
	float *v1,
	float *v2,
	float *v3,
	float *v4,
	float *v5,
	float *v6,
	float *v7
	) {
	ERR_INIT

	const size_t itensitiesInterleavedSize = maxVecSize * AVX2_FLOAT_REGISTER_SIZE;
	const size_t itensitiesInterleavedSizeAlignas = calculateNextAlignedBlockSize(
		itensitiesInterleavedSize,
		AVX2_ALIGNAS_SIZE
		);
	alignas(AVX2_ALIGNAS_SIZE) float itensitiesInterleaved[itensitiesInterleavedSizeAlignas] = {0};

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
		 itensitiesInterleaved
		);

	__m256 registerMaxes = _mm256_setzero_ps();
	__m256 inPeakMask = _mm256_setzero_ps();
	for (int i = 0; i < itensitiesInterleavedSizeAlignas; i++) {

		const __m256 frameIndexRegisterCurrent
			= _mm256_load_ps(itensitiesInterleaved + (i * AVX2_FLOAT_REGISTER_SIZE));

		__m256 registerMaxesComparison = _mm256_cmp_ps(frameIndexRegisterCurrent, registerMaxes, _CMP_GT_OQ);


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

void AVXUtils::cosineSimilarityAVXParallel(
	const float *arrRef,
	const float *array0,
	const float *array1,
	const float *array2,
	const float *array3,
	const float *array4,
	const float *array5,
	const float *array6,
	const float *array7,
	size_t length,
	float* cosineSimResultsAligned
	) {

	__m256 dotProductsBuilder = _mm256_setzero_ps();
	__m256 magnitudesBuilder = _mm256_setzero_ps();
	__m256 magnitudesBuilderRef = _mm256_setzero_ps();

	for (int i = 0; i < length; i ++) {

		const __m256 ref = _mm256_set1_ps(arrRef[i]);
		const __m256 magnitudeRef = _mm256_mul_ps(ref, ref);
		magnitudesBuilderRef = _mm256_add_ps(magnitudesBuilderRef, magnitudeRef);

		const __m256 vi = _mm256_set_ps(
			array7[i], array6[i], array5[i], array4[i], array3[i], array2[i], array1[i], array0[i]
			);

		const __m256 producti = _mm256_mul_ps(ref, vi);
		dotProductsBuilder = _mm256_add_ps(dotProductsBuilder, producti);

		const __m256 magnitudei = _mm256_mul_ps(vi, vi);
		magnitudesBuilder = _mm256_add_ps(magnitudesBuilder, magnitudei);
	}

	const __m256 magnitudesSqrt = _mm256_sqrt_ps(magnitudesBuilder);
	const __m256 magnitudeRefSqrt = _mm256_sqrt_ps(magnitudesBuilderRef);

	const __m256 denomonator = _mm256_mul_ps(magnitudesSqrt, magnitudeRefSqrt);

	__m256 resultVec = _mm256_div_ps(dotProductsBuilder, denomonator);

	const __m256 nanMask = _mm256_cmp_ps(resultVec, resultVec, _CMP_UNORD_Q);
	const __m256 zero_vec = _mm256_setzero_ps();
	resultVec = _mm256_blendv_ps(resultVec, zero_vec, nanMask);

	_mm256_store_ps(cosineSimResultsAligned, resultVec);
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








