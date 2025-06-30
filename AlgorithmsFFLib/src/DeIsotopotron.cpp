//
// Created by andrewnichols on 6/29/25.
//

#include "DeIsotopotron.h"

#include "MathUtils.h"
#include "ParallelUtils.h"

#include <cpuid.h>
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

	void copyFloatAVX512(
		const float* src,
		float* dst,
		size_t count
		){
			constexpr size_t simdWidth = 16;
			for (size_t i = 0; i + simdWidth < count; i += simdWidth) {
				__m512 v = _mm512_loadu_ps(src + i);
				_mm512_store_ps(dst + i, v);
			}
	}

	void copyFloatAVX2(
		const float* src,
		float* dst,
		size_t count
		) {
		constexpr size_t simdWidth = 8;
		for (size_t i = 0; i + simdWidth <= count; i += simdWidth) {
			__m256 v = _mm256_loadu_ps(src + i);
			_mm256_storeu_ps(dst + i, v);
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
		const __m256& avx256Mz,
		float mzMin,
		float mzMax,
		float subtractVal,
		__m256* avx256Intensity
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

}//namespace
Err DeIsotopotron::deisotopeTandemScan(
	float extractionTolerancePPM,
	ScanPoints *scanPoints
	) {

	ERR_INIT

	QPair<QVector<float>, QVector<float>> scanPointsSplit = ParallelUtils::unZipPointFF(*scanPoints);
	const QVector<float> &mzVals = scanPointsSplit.first;
	const QVector<float> &intensityVals = scanPointsSplit.second;

	if (supportsAVX512()) {

		constexpr int avx512FloatSize = 16;
		const int vecSizeBuffered = calculateBufferedVecSizeAVX(scanPoints->size(), avx512FloatSize);
		const int avx512RegistersCount = vecSizeBuffered / avx512FloatSize;

		alignas(64) float mzValsRaw[vecSizeBuffered] = {0.0f};
		copyFloatAVX512(mzVals.data(), mzValsRaw, mzVals.size());

		alignas(64) float intensityValsRaw[vecSizeBuffered] = {0.0f};
		copyFloatAVX512(intensityVals.data(), intensityValsRaw, intensityVals.size());

		__m512 mzRegisterCurrent = _mm512_setzero_ps();
		__m512 mzRegisterNext = _mm512_setzero_ps();
		__m512 intensityRegisterCurrent = _mm512_setzero_ps();
		__m512 intensityRegisterNext = _mm512_setzero_ps();
		printAVX512Float(intensityRegisterCurrent);

		int avxIndexCurrent = -1;
		for (int i = 0; i < mzVals.size(); i++) {

			const int avxIndex = i / avx512FloatSize;
			if (avxIndex != avxIndexCurrent) {

				if (i > 0) {

					// printAVX512Float(intensityRegisterCurrent);

					// _mm512_storeu_ps(
					// 	&intensityValsRaw + (avxIndexCurrent * avx512FloatSize),
					// 	intensityRegisterCurrent
					// 	);

					// _mm512_store_ps(
					// 	intensityValsRaw + ((avxIndexCurrent + 1) * avx512FloatSize),
					// 	intensityRegisterNext
					// 	);
				}


				mzRegisterCurrent = _mm512_load_ps(mzValsRaw + (avxIndex * avx512FloatSize));
				intensityRegisterCurrent = _mm512_load_ps(intensityValsRaw + (avxIndex * avx512FloatSize));

				if (avxIndex + 1 < avx512RegistersCount) {
					mzRegisterNext = _mm512_load_ps(mzValsRaw + ((avxIndex + 1) * avx512FloatSize));
					intensityRegisterNext = _mm512_load_ps(intensityValsRaw + ((avxIndex + 1) * avx512FloatSize));
				}

				avxIndexCurrent = avxIndex;
			}

			const float mzValCurrent = mzValsRaw[i];
			const float mzValCurrentCharge1 = mzValCurrent + static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF);
			const float mzValCurrentCharge2 = mzValCurrent + static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF) / 2.0f;
			const float mzTol = MathUtils::calculatePPM(mzValCurrentCharge1, extractionTolerancePPM);
			const float intensityValCurrent = intensityValsRaw[i];

			if (MathUtils::tZero(intensityValCurrent)) {
				continue;
			}

			subtractValInRangeFromIntensityAtGivenMz(
				mzRegisterCurrent,
				mzValCurrentCharge1 - mzTol,
				mzValCurrentCharge1 + mzTol,
				intensityValCurrent,
				&intensityRegisterCurrent
				);

			subtractValInRangeFromIntensityAtGivenMz(
				mzRegisterCurrent,
				mzValCurrentCharge2 - mzTol,
				mzValCurrentCharge2 + mzTol,
				intensityValCurrent,
				&intensityRegisterCurrent
				);

			if (avxIndexCurrent + 1 < avx512RegistersCount) {

				subtractValInRangeFromIntensityAtGivenMz(
					mzRegisterNext,
					mzValCurrentCharge1 - mzTol,
					mzValCurrentCharge1 + mzTol,
					intensityValCurrent,
					&intensityRegisterNext
					);

				subtractValInRangeFromIntensityAtGivenMz(
					mzRegisterNext,
					mzValCurrentCharge2 - mzTol,
					mzValCurrentCharge2 + mzTol,
					intensityValCurrent,
					&intensityRegisterNext
					);
			}

		}
	}

	else {
		constexpr int avx256FloatSize = 8;
		const int vecSizeBuffered = calculateBufferedVecSizeAVX(scanPoints->size(), avx256FloatSize);
		const int avx256RegistersCount = vecSizeBuffered / avx256FloatSize;

		alignas(32) float mzValsRaw[vecSizeBuffered] = {0.0f};
		copyFloatAVX2(mzVals.data(), mzValsRaw, mzVals.size());

		alignas(32) float intensityValsRaw[vecSizeBuffered] = {0.0f};
		copyFloatAVX2(intensityVals.data(), intensityValsRaw, intensityVals.size());

		__m256 mzRegisterCurrent = _mm256_setzero_ps();
		__m256 mzRegisterNext = _mm256_setzero_ps();
		__m256 intensityRegisterCurrent = _mm256_setzero_ps();
		__m256 intensityRegisterNext = _mm256_setzero_ps();

		int avxIndexCurrent = -1;
		for (int i = 0; i < mzVals.size(); i++) {
		    const int avxIndex = i / avx256FloatSize;
		    if (avxIndex != avxIndexCurrent) {
		        if (i > 0) {

		            _mm256_store_ps(
		                &intensityValsRaw[(avxIndexCurrent * avx256FloatSize)],
		                intensityRegisterCurrent
		            );

		        	if (avxIndexCurrent + 1 < avx256RegistersCount) {

remove the following and use another for (int i = 0; i < mzVals.size(); i++) with intensityValsRaw[((avxIndex + 1) * avx256FloatSize)],

		        		// _mm256_store_ps(
		        		//     &intensityValsRaw[((avxIndex + 1) * avx256FloatSize)],
		        		//     intensityRegisterNext
		        		// );
		        	}


		        }

		    	avxIndexCurrent = avxIndex;

		        mzRegisterCurrent = _mm256_load_ps(mzValsRaw + (avxIndexCurrent * avx256FloatSize));
		        intensityRegisterCurrent = _mm256_load_ps(intensityValsRaw + (avxIndexCurrent * avx256FloatSize));

		        if (avxIndex + 1 < avx256RegistersCount) {
		            mzRegisterNext = _mm256_load_ps(mzValsRaw + ((avxIndexCurrent + 1) * avx256FloatSize));
		            intensityRegisterNext = _mm256_load_ps(intensityValsRaw + ((avxIndexCurrent + 1) * avx256FloatSize));
		        }

		    }

		    const float mzValCurrent = mzValsRaw[i];
		    const float mzValCurrentIsotopeCharge1 = mzValCurrent + (static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF) / 1.0f);
		    const float mzValCurrentIsotopeCharge2 = mzValCurrent + (static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF) / 2.0f);
		    const float mzTol = MathUtils::calculatePPM(mzValCurrent, extractionTolerancePPM);
		    const float intensityValCurrent = intensityValsRaw[i];

		    if (MathUtils::tZero(intensityValCurrent)) {
		        continue;
		    }

		    subtractValInRangeFromIntensityAtGivenMzAVX2(
		        mzRegisterCurrent,
		        mzValCurrentIsotopeCharge1 - mzTol,
		        mzValCurrentIsotopeCharge1 + mzTol,
		        intensityValCurrent,
		        &intensityRegisterCurrent
		    );

		    subtractValInRangeFromIntensityAtGivenMzAVX2(
		        mzRegisterCurrent,
		        mzValCurrentIsotopeCharge2 - mzTol,
		        mzValCurrentIsotopeCharge2 + mzTol,
		        intensityValCurrent,
		        &intensityRegisterCurrent
		    );

		    if (avxIndexCurrent + 1 < avx256RegistersCount) {
		        subtractValInRangeFromIntensityAtGivenMzAVX2(
		            mzRegisterNext,
		            mzValCurrentIsotopeCharge1 - mzTol,
		            mzValCurrentIsotopeCharge1 + mzTol,
		            intensityValCurrent,
		            &intensityRegisterNext
		        );

		        subtractValInRangeFromIntensityAtGivenMzAVX2(
		            mzRegisterNext,
		            mzValCurrentIsotopeCharge2 - mzTol,
		            mzValCurrentIsotopeCharge2 + mzTol,
		            intensityValCurrent,
		            &intensityRegisterNext
		        );
		    }
		}


		for (int i = 0; i < mzVals.size(); i++) {
			qDebug() << mzVals[i] << mzValsRaw[i] << intensityVals[i] << intensityValsRaw[i];
		}




	}






	ERR_RETURN
}
