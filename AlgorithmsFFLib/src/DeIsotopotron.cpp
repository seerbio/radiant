//
// Created by andrewnichols on 6/29/25.
//

#include "DeIsotopotron.h"

#include "MathUtils.h"
#include "ParallelUtils.h"

#include <cpuid.h>
#include <cstring>
#include <immintrin.h>


namespace {

	int calculateBufferedVecSizeAVX(
		int vecSizeOG,
		int desiredDivisor
		) {

		while (vecSizeOG % desiredDivisor != 0) {
			vecSizeOG++;
		}

		return vecSizeOG;
	}


	void copyFloatAVX2(
		const float* src,
		float* dst,
		size_t count
		) {
		constexpr size_t simdWidth = 8;
		size_t i = 0;

		for (; i + simdWidth <= count; i += simdWidth) {
			__m256 v = _mm256_loadu_ps(src + i);
			_mm256_storeu_ps(dst + i, v);
		}

		if (i < count) {
			std::memcpy(dst + i, src + i, (count - i) * sizeof(float));
		}

	}

	void subtractValInRangeFromIntensityAtGivenMz(
		const __m512 &avx512Mz,
		float mzMin,
		float mzMax,
		float subtractVal,
		__m512 *avx512Intensity
		) {

		const __m512 mzValMin = _mm512_set1_ps(mzMin);
		const __m512 mzValMax = _mm512_set1_ps(mzMax);
		const __m512 subtractVals = _mm512_set1_ps(subtractVal);

		__mmask16 maskMin = _mm512_cmp_ps_mask(avx512Mz, mzValMin, _CMP_GE_OQ);
		__mmask16 maskMax = _mm512_cmp_ps_mask(avx512Mz, mzValMax, _CMP_LE_OQ);
		__mmask16 maskCombined = maskMin & maskMax;
		*avx512Intensity = _mm512_mask_sub_ps(
			*avx512Intensity,
			maskCombined,
			*avx512Intensity,
			subtractVals
			);
	}

	void subtractValInRangeFromIntensityAtGivenMzAVX2(
		const __m256 &avx256Mz,
		float mzMin,
		float mzMax,
		float subtractVal,
		__m256 *avx256Intensity
		) {
			const __m256 mzValMin = _mm256_set1_ps(mzMin);
			const __m256 mzValMax = _mm256_set1_ps(mzMax);
			const __m256 subtractVals = _mm256_set1_ps(subtractVal);

			__m256 maskMin = _mm256_cmp_ps(avx256Mz, mzValMin, _CMP_GE_OQ); // Get comparisons
			__m256 maskMax = _mm256_cmp_ps(avx256Mz, mzValMax, _CMP_LE_OQ); // Check if <= mzMax

			__m256 maskCombined = _mm256_and_ps(maskMin, maskMax);

			*avx256Intensity = _mm256_blendv_ps(
				*avx256Intensity,
				_mm256_sub_ps(*avx256Intensity, subtractVals),
				maskCombined
			);
		}

	static void printAVX(const __m256 &avx) {
		alignas(64) float values[8];
		_mm256_store_ps(values, avx);
		QDebug debug = qDebug().nospace();
		for (int j = 0; j < 8; j++) {
			debug << MathUtils::pRound(values[j], 3) << ", ";
		}
	}

	void printAVX512Float(const __m512 &avx) {
		alignas(64) float values[16];
		_mm512_storeu_ps(values, avx);
		for (int j = 0; j < 16; j++) {
			std::cout << values[j] << ", ";     // Print each value
		}
		std::cout << std::endl;
	}

	void printAVX2Float(const __m256& avx) {
		alignas(32) float values[8];   // Each AVX2 register holds 8 floats
		_mm256_storeu_ps(values, avx); // Store the contents of the AVX2 register into an array
		for (int j = 0; j < 8; j++) {
			std::cout << values[j] << ", ";     // Print each value
		}
		std::cout << std::endl;
	}

	bool supportsAVX512() {
		unsigned int eax, ebx, ecx, edx;
		__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx);
		return (ebx & (1 << 16)) != 0; // AVX512F bit in EBX
	}

	float horizontalMin(__m256 vec) {
		// Extract the 8 floats into a 256-bit vector
		__m128 low = _mm256_castps256_ps128(vec); // Lower 4 floats
		__m128 high = _mm256_extractf128_ps(vec, 1); // Upper 4 floats
		__m128 min4 = _mm_min_ps(low, high); // Element-wise minimum of low and high

		// Reduce the result further
		__m128 min2 = _mm_min_ps(min4, _mm_movehl_ps(min4, min4)); // Compare first two and last two floats
		__m128 min1 = _mm_min_ps(min2, _mm_shuffle_ps(min2, min2, 0x1)); // Compare adjacent floats

		return _mm_cvtss_f32(min1); // Extract the smallest float
	}

	float horizontalMax(__m256 vec) {
		__m128 low = _mm256_castps256_ps128(vec); // Lower 4 floats
		__m128 high = _mm256_extractf128_ps(vec, 1); // Upper 4 floats
		__m128 max4 = _mm_max_ps(low, high); // Element-wise maximum of low and high

		__m128 max2 = _mm_max_ps(max4, _mm_movehl_ps(max4, max4)); // Compare first two and last two floats
		__m128 max1 = _mm_max_ps(max2, _mm_shuffle_ps(max2, max2, 0x1)); // Compare adjacent floats

		return _mm_cvtss_f32(max1); // Extract the largest float
	}

	__m256 fillNegativesWithZeros(__m256 input) {
		__m256 zeros = _mm256_setzero_ps();
		__m256 mask = _mm256_cmp_ps(input, zeros, _CMP_LT_OS);
		__m256 result = _mm256_blendv_ps(input, zeros, mask);
		return result;
	}

}//namespace
Err DeIsotopotron::deisotopeTandemScan(
	float extractionTolerancePPM,
	ScanPoints *scanPoints
	) {

	ERR_INIT

	QPair<QVector<float>, QVector<float>> scanPointsSplit = ParallelUtils::unZipPointFF(*scanPoints);
	const QVector<float> &mzVals = scanPointsSplit.first;
	const QVector<float> &intensityVals = scanPointsSplit.second;

	constexpr int avx256FloatSize = 8;
	const int vecSizeBuffered = calculateBufferedVecSizeAVX(scanPoints->size(), avx256FloatSize);

	alignas(32) float mzValsRaw[vecSizeBuffered] = {0.0f};
	copyFloatAVX2(mzVals.data(), mzValsRaw, mzVals.size());

	alignas(32) float intensityValsRaw[vecSizeBuffered] = {0.0f};
	copyFloatAVX2(intensityVals.data(), intensityValsRaw, intensityVals.size());

	float mzRegisterMin = 0.0f;
	int currentRegisterIndex = 0;
	__m256 mzRegisterCurrent = _mm256_load_ps(mzValsRaw + currentRegisterIndex);
	__m256 intensityzRegisterCurrent = _mm256_load_ps(intensityValsRaw + currentRegisterIndex);
	float mzRegisterMax = horizontalMax(mzRegisterCurrent);
	for (int i = 0; i < mzVals.size(); i++) {

		const float mzVal = mzValsRaw[i];

		float mzValIso1 = mzVal + static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF * 0.5);
		float mzValIso2 = mzVal + static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF);
		const float mzTol = MathUtils::calculatePPM(mzVal, extractionTolerancePPM);

		if (
			(mzValIso1 + mzTol < mzRegisterMin && mzValIso2 + mzTol < mzRegisterMin) ||
			MathUtils::tSame(mzVal, mzRegisterMin) || MathUtils::tSame(mzVal, mzRegisterMax)
			) {
			continue;
		}

		while (mzValIso1 - mzTol > mzRegisterMax && mzValIso2 - mzTol > mzRegisterMax) {
			_mm256_store_ps(intensityValsRaw + currentRegisterIndex, intensityzRegisterCurrent);
			currentRegisterIndex += avx256FloatSize;
			mzRegisterCurrent = _mm256_load_ps(mzValsRaw + currentRegisterIndex);
			intensityzRegisterCurrent = _mm256_load_ps(intensityValsRaw + currentRegisterIndex);
			mzRegisterMin = mzRegisterMax;
			mzRegisterMax = horizontalMax(mzRegisterCurrent);
		}

		subtractValInRangeFromIntensityAtGivenMzAVX2(
			mzRegisterCurrent,
			mzValIso1 - mzTol,
			mzValIso1 + mzTol,
			intensityValsRaw[i],
			&intensityzRegisterCurrent
			);

		subtractValInRangeFromIntensityAtGivenMzAVX2(
			mzRegisterCurrent,
			mzValIso2 - mzTol,
			mzValIso2 + mzTol,
			intensityValsRaw[i],
			&intensityzRegisterCurrent
			);

		// qDebug() << mzValIso1 - mzTol << mzValIso1 + mzTol;
		// qDebug() << mzValIso2 - mzTol << mzValIso2 + mzTol;
		// printAVX(mzRegisterCurrent);
		// printAVX(intensityzRegisterCurrent);
		// intensityzRegisterCurrent = fillNegativesWithZeros(intensityzRegisterCurrent);
	}

	scanPoints->clear();
	for (int i = 0; i < mzVals.size(); i++) {
		const float intensityVal = intensityValsRaw[i];
		if (intensityVal < 0) {
			continue;
		}
		scanPoints->push_back({mzValsRaw[i], intensityValsRaw[i]});
	}

	ERR_RETURN
}
