#include "AVXUtils.h"


#include <QtTest>


class AVXUtilsTests : public QObject
{
    Q_OBJECT

public:
    AVXUtilsTests() = default;
    ~AVXUtilsTests() override = default;

private slots:
    static void copyAVXTest();
	static void calculateNextAlignedBlockSizeTest();
	static void printAVXFloatTest();
	static void printAVXMaskTest();
	static void convolveWithKernelAVXFloatTest();
	static void splitAVXUInt16to32Test();
	static void subtractArraysAVX2Test();
	static void isNByteAlignedTest();
	static void replaceArrayValuesAVXTest();
	static void cosineSimilarityAVXTest();
	static void cosineSimilarityAVXTest2();
	static void maxFloatTest();
	static void isAllOnesTest();
	static void interleaveVectorsTest();
	static void separateInterleavedVectorsTest();
	static void findApexesEightVecsTest();
	static void log256Test();
	static void cosineSimilarityAVXParallelTest();
	static void sumParallelTest();
	static void meanParallelTest();
	static void stDevParallelTest();
};


void AVXUtilsTests::copyAVXTest() {

	ERR_INIT

	size_t size = 96;
	alignas(32) float src[size] = {0};
	alignas(32) float dst[size] = {0};

	for (size_t i = 0; i < size; i++) {
		src[i] = static_cast<float>(i);
	}

	e = AVXUtils::copyAVXFloatToAligned(src, dst, size);
	QCOMPARE(e, Err::eNoError);

	for (size_t i = 0; i < size; i++) {
		QCOMPARE(static_cast<int>(dst[i]), static_cast<int>(src[i]));
	}

	alignas(32) float src1[size - 2] = {0};
	alignas(32) float dst1[size] = {0};

	for (size_t i = 0; i < size - 2; i++) {
		src1[i] = static_cast<float>(i);
	}
	e = AVXUtils::copyAVXFloatToAligned(src1, dst1, size-2);
	QCOMPARE(e, Err::eNoError);

	for (size_t i = 0; i < size - 2; i++) {
		QCOMPARE(static_cast<int>(dst[i]), static_cast<int>(src1[i]));
	}
	QCOMPARE(static_cast<int>(dst1[size - 3]), 93);
	QCOMPARE(static_cast<int>(dst1[size - 2]), 0);
	QCOMPARE(static_cast<int>(dst1[size - 1]), 0);
}

void AVXUtilsTests::calculateNextAlignedBlockSizeTest() {

	const size_t nextAlignedSize
		= AVXUtils::calculateNextAlignedBlockSize(666, AVXUtils::AVX2_FLOAT_REGISTER_SIZE);

	QCOMPARE(nextAlignedSize, 672);

}

void AVXUtilsTests::printAVXFloatTest() {
	const std::vector<float> v = {1,2,3,4,5,6,7,8};
	const __m256 data = _mm256_loadu_ps(v.data());
	AVXUtils::printAVXFloat(data);
}

void AVXUtilsTests::printAVXMaskTest() {

	const __m256 mzValMin = _mm256_set1_ps(4.0f);

	const std::vector<float> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	const __m256 data = _mm256_loadu_ps(v.data());
	const __m256 mask = _mm256_cmp_ps(data, mzValMin, _CMP_GE_OQ);

	AVXUtils::printMask(mask, 8);
}

void AVXUtilsTests::convolveWithKernelAVXFloatTest() {

	ERR_INIT

	constexpr size_t size = 21;

	QVector<float> v0(size, 0);
	QVector<float> v1(size, 1);
	QVector<float> v2(size, 2);
	QVector<float> v3(size, 3);
	QVector<float> v4(size, 4);
	QVector<float> v5(size, 5);
	QVector<float> v6(size, 6);
	QVector<float> v7(size, 7);

	QVector<float> vAsc(size);
	std::iota(vAsc.begin(), vAsc.end(), 1.0f);

	QVector<float> oneHotMid(size, 0);
	oneHotMid[oneHotMid.size() / 2 - 1] = 1;

	QVector<float> oneHotMidCopy(size, 0);
	oneHotMidCopy[oneHotMidCopy.size() / 2 - 1] = 1;


	QVector<float> oneHotLast(size, 0);
	oneHotLast[oneHotLast.size() - 1] = 1;

	// const QVector<float> kernel = {0.25, 0.5, 0.25};
	// const QVector<float> kernel = {0.1, 0.2, 0.4, 0.2, 0.1};
	const QVector<float> kernel = {0.1, 1, 2, 4, 2, 1, 0.1};

	e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
		kernel,
		size,
		oneHotLast.data(),
		v0.data(),
		v0.data(),
		v0.data(),
		oneHotMid.data(),
		v0.data(),
		v0.data(),
		v0.data()
		);

		QCOMPARE(static_cast<int>(oneHotMid[oneHotMid.size() / 2 - 1]), 4);
}

void AVXUtilsTests::splitAVXUInt16to32Test() {

	ERR_INIT

	alignas(32) uint16_t inputValues[16] = {
		1, 2, 3, 4, 5, 6, 7, 8,
		9, 10, 11, 12, 13, 14, 15, 16
	};

	__m256i input = _mm256_load_si256(reinterpret_cast<const __m256i*>(inputValues));

	__m256i int32Output1;
	__m256i int32Output2;
	e = AVXUtils::splitAVXUInt16to32(
		input,
		&int32Output1,
		&int32Output2
		);

	// AVXUtils::printAVXInt32(int32Output1);
	// AVXUtils::printAVXInt32(int32Output2);

	__m256i int32OutputRecombine = _mm256_packs_epi32(int32Output1, int32Output2);
	// AVXUtils::printAVXInt16(int32OutputRecombine);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) int16_t values[AVXUtils::AVX2_INT16_REGISTER_SIZE];
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(values), int32OutputRecombine);

	for (int i = 0; i < AVXUtils::AVX2_INT16_REGISTER_SIZE; i++) {
		QCOMPARE(values[i], inputValues[i]);
	}

}

void AVXUtilsTests::subtractArraysAVX2Test() {

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float arr1[16];
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float arr2[16];
	for (int i = 0; i < 16; i++) {
		arr1[i] = static_cast<float>(i);
		arr2[i] = static_cast<float>(i);
	}
	AVXUtils::subtractArraysAVX2(
		arr1,
		arr2,
		16,
		false
		);

	for (int i = 0; i < 16; i++) {
		QVERIFY(MathUtils::tZero(arr1[i]));
	}

}

void AVXUtilsTests::isNByteAlignedTest() {

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float unaligned[16];
	float* unalignedPtr = reinterpret_cast<float*>(
		reinterpret_cast<uint8_t*>(unaligned) + 4
	);
	const bool isAligned = AVXUtils::isNByteAligned(unalignedPtr, AVXUtils::AVX2_ALIGNAS_SIZE);
	QCOMPARE(isAligned, false);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float aligned[16];
	const bool isAligned2 = AVXUtils::isNByteAligned(aligned, AVXUtils::AVX2_ALIGNAS_SIZE);
	QCOMPARE(isAligned2, true);

}

void AVXUtilsTests::replaceArrayValuesAVXTest() {

	alignas(32) float arr1[8] = {0, 100, 200, 0, 0, 100, 200, 0};
	__m256 values = _mm256_load_ps(arr1);

	AVXUtils::replaceArrayValuesAVXGreaterThan(0.1, 1.0, values);

	_mm256_store_ps(arr1, values);

	alignas(32) float expectedResult[8] = {0, 1, 1, 0, 0, 1, 1, 0};
	for (int i = 0; i < 8; i++) {
		QCOMPARE(static_cast<int>(arr1[i]), static_cast<int>(expectedResult[i]));
	}
}

void AVXUtilsTests::cosineSimilarityAVXTest() {

	const __m256 a = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
	const __m256 b = _mm256_set_ps(0, 1, 2, 3, 4, 5, 6, 7);
	const __m256 c = _mm256_set_ps(1, 1, 1, 1, 0, 0, 0, 0);
	const __m256 d = _mm256_set_ps(2, 0, 0, 0, 1, 1, 0, 0);

	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(a, b), 0.4f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(a, a), 1.0f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(_mm256_setzero_ps(), a), 0.0f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(c, d), 0.408248f));

}

void AVXUtilsTests::cosineSimilarityAVXTest2() {

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float z[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float a[8] = {7, 6, 5, 4, 3, 2, 1, 0};
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float b[8] = {0, 1, 2, 3, 4, 5, 6, 7};
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float c[8] = {1, 1, 1, 1, 0, 0, 0, 0};
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float d[8] = {2, 0, 0, 0, 1, 1, 0, 0};

	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(a, b, 8), 0.4f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(a, a, 8), 1.0f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(z, a, 8), 0.0f));
	QVERIFY(MathUtils::tSame(AVXUtils::cosineSimilarityAVX(c, d, 8), 0.408248f));

}

void AVXUtilsTests::maxFloatTest() {
	__m256 reg = _mm256_set_ps(1.0, -2.0, 22.0, 11.0, 133.0, 3.0, 666.0, 665.9);
	const float regMax = AVXUtils::maxFloat(reg);
	QCOMPARE(static_cast<int>(regMax), 666);
}

void AVXUtilsTests::isAllOnesTest() {

	const __m256 six = _mm256_set1_ps(6);
	const __m256 eight = _mm256_set1_ps(8);
	const __m256 mix = _mm256_set_ps(2.0, 7.0, 2.0, 7.0, 2.0, 7.0, 2.0, 7.0);

	const __m256 maskPos = _mm256_cmp_ps(mix, eight, _CMP_LT_OQ);
	const __m256 maskNeg = _mm256_cmp_ps(mix, six, _CMP_LT_OQ);

	const bool posTest = AVXUtils::isAllOnes(maskPos);
	QCOMPARE(posTest, true);

	const bool negTest = AVXUtils::isAllOnes(maskNeg);
	QCOMPARE(negTest, false);
}

void AVXUtilsTests::interleaveVectorsTest() {

	QVector<float> ones(8, 1);
	QVector<float> twos(8, 2);
	QVector<float> threes(8, 3);
	QVector<float> fours(8, 4);
	QVector<float> fives(8, 5);
	QVector<float> sixes(8, 6);
	QVector<float> sevens(8, 7);
	QVector<float> eights(8, 8);

	int padding = 0;
	int resultSize = AVXUtils::AVX2_FLOAT_REGISTER_SIZE * (ones.size() + (padding * 2));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float result[resultSize] = {0};
	AVXUtils::interleaveVectors(
		ones.size(),
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data(),
		result
		);

	for (int i = 0; i < ones.size(); i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE) {
		QCOMPARE(static_cast<int>(result[i]), 1);
		QCOMPARE(static_cast<int>(result[i+1]), 2);
		QCOMPARE(static_cast<int>(result[i+2]), 3);
		QCOMPARE(static_cast<int>(result[i+3]), 4);
		QCOMPARE(static_cast<int>(result[i+4]), 5);
		QCOMPARE(static_cast<int>(result[i+5]), 6);
		QCOMPARE(static_cast<int>(result[i+6]), 7);
		QCOMPARE(static_cast<int>(result[i+7]), 8);
	}

	padding = 2;
	resultSize = AVXUtils::AVX2_FLOAT_REGISTER_SIZE * (ones.size() + (padding * 2));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultPad[resultSize] = {0};
	AVXUtils::interleaveVectors(
		ones.size(),
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data(),
		resultPad
		);

	for (int i = 0; i < resultSize; i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE) {

		if (i < (2 * AVXUtils::AVX2_FLOAT_REGISTER_SIZE) || i >= resultSize - (2 * AVXUtils::AVX2_FLOAT_REGISTER_SIZE)) {
			for (int j = 0; j < 8; j++) {
				QCOMPARE(static_cast<int>(resultPad[i + j]), 0);
			}
			continue;
		}
		QCOMPARE(static_cast<int>(resultPad[i]), 1);
		QCOMPARE(static_cast<int>(resultPad[i+1]), 2);
		QCOMPARE(static_cast<int>(resultPad[i+2]), 3);
		QCOMPARE(static_cast<int>(resultPad[i+3]), 4);
		QCOMPARE(static_cast<int>(resultPad[i+4]), 5);
		QCOMPARE(static_cast<int>(resultPad[i+5]), 6);
		QCOMPARE(static_cast<int>(resultPad[i+6]), 7);
		QCOMPARE(static_cast<int>(resultPad[i+7]), 8);
	}
}

void AVXUtilsTests::separateInterleavedVectorsTest() {

	QVector<float> ones(8, 1);
	QVector<float> twos(8, 2);
	QVector<float> threes(8, 3);
	QVector<float> fours(8, 4);
	QVector<float> fives(8, 5);
	QVector<float> sixes(8, 6);
	QVector<float> sevens(8, 7);
	QVector<float> eights(8, 8);

	int padding = 0;
	int resultSize = AVXUtils::AVX2_FLOAT_REGISTER_SIZE * (ones.size() + (padding * 2));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float result[resultSize] = {0};
	AVXUtils::interleaveVectors(
		ones.size(),
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data(),
		result
		);

	AVXUtils::separateInterleavedVectors(
		result,
		resultSize,
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data()
		);

	for (int i = 0; i < 8; i ++) {
		QCOMPARE(static_cast<int>(ones[i]), 1);
		QCOMPARE(static_cast<int>(twos[i]), 2);
		QCOMPARE(static_cast<int>(threes[i]), 3);
		QCOMPARE(static_cast<int>(fours[i]), 4);
		QCOMPARE(static_cast<int>(fives[i]), 5);
		QCOMPARE(static_cast<int>(sixes[i]), 6);
		QCOMPARE(static_cast<int>(sevens[i]), 7);
		QCOMPARE(static_cast<int>(eights[i]), 8);
	}

	// qDebug()
	// << resultPad[i]
	// << resultPad[i+1]
	// << resultPad[i+2]
	// << resultPad[i+3]
	// << resultPad[i+4]
	// << resultPad[i+5]
	// << resultPad[i+6]
	// << resultPad[i+7]
	// ;

	padding = 3;
	resultSize = AVXUtils::AVX2_FLOAT_REGISTER_SIZE * (ones.size() + (padding * 2));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultPad[resultSize] = {0};
	AVXUtils::interleaveVectors(
		ones.size(),
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data(),
		resultPad
		);

	AVXUtils::separateInterleavedVectors(
		resultPad,
		resultSize,
		padding,
		ones.data(),
		twos.data(),
		threes.data(),
		fours.data(),
		fives.data(),
		sixes.data(),
		sevens.data(),
		eights.data()
		);

	for (int i = 0; i < 8; i ++) {
		QCOMPARE(static_cast<int>(ones[i]), 1);
		QCOMPARE(static_cast<int>(twos[i]), 2);
		QCOMPARE(static_cast<int>(threes[i]), 3);
		QCOMPARE(static_cast<int>(fours[i]), 4);
		QCOMPARE(static_cast<int>(fives[i]), 5);
		QCOMPARE(static_cast<int>(sixes[i]), 6);
		QCOMPARE(static_cast<int>(sevens[i]), 7);
		QCOMPARE(static_cast<int>(eights[i]), 8);
	}
}

void AVXUtilsTests::findApexesEightVecsTest() {

	const int size = 8;

	QVector<float> v0(size, 0);
	QVector<float> v1(size, 0);
	QVector<float> v2(size, 0);
	QVector<float> v3(size, 0);
	QVector<float> v4(size, 0);
	QVector<float> v5(size, 0);
	QVector<float> v6(size, 0);
	QVector<float> v7(size, 0);
	QVector<float> vApexes0(size, 0);
	QVector<float> vApexes1(size, 0);
	QVector<float> vApexes2(size, 0);
	QVector<float> vApexes3(size, 0);
	QVector<float> vApexes4(size, 0);
	QVector<float> vApexes5(size, 0);
	QVector<float> vApexes6(size, 0);
	QVector<float> vApexes7(size, 0);

	v0[0] = 0;
	v1[1] = 1;
	v2[2] = 2;
	v3[3] = 3;
	v4[4] = 4;
	v5[5] = 5;
	v6[6] = 6;
	v7[7] = 7;

	v0[0] = 0;
	v1[0] = 1;
	v2[1] = 1;
	v3[2] = 1;
	v4[3] = 1;
	v5[4] = 1;
	v6[5] = 1;
	v7[6] = 1;

	const size_t masterVecApexesSize = AVXUtils::calculateNextAlignedBlockSize(
		size * AVXUtils::AVX2_FLOAT_REGISTER_SIZE,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float masterVecApexes[masterVecApexesSize];

	ERR_INIT
	e = AVXUtils::findApexesEightVecs(
		size,
		v0.data(),
		v1.data(),
		v2.data(),
		v3.data(),
		v4.data(),
		v5.data(),
		v6.data(),
		v7.data(),
		masterVecApexes
		);

	AVXUtils::separateInterleavedVectors(
		masterVecApexes,
		masterVecApexesSize,
		0,
		vApexes0.data(),
		vApexes1.data(),
		vApexes2.data(),
		vApexes3.data(),
		vApexes4.data(),
		vApexes5.data(),
		vApexes6.data(),
		vApexes7.data()
		);

	QCOMPARE(static_cast<int>(v2[2]), 2);
	QCOMPARE(static_cast<int>(v3[3]), 3);
	QCOMPARE(static_cast<int>(v4[4]), 4);
	QCOMPARE(static_cast<int>(v5[5]), 5);
	QCOMPARE(static_cast<int>(v6[6]), 6);
}

void AVXUtilsTests::log256Test() {

	__m256 v0 = _mm256_set_ps(2, 4, 8, 16, 32, 64, 128, 256);

	__m256 logLn = AVXUtils::log256(v0);

	float v[8];
	_mm256_store_ps(v, logLn);

	int e[8] = {5, 4, 3, 3, 2, 1, 1, 0};

	for (int i = 0; i < 8; i ++) {
		QCOMPARE(static_cast<int>(v[i]), e[i]);
	}
}

void AVXUtilsTests::cosineSimilarityAVXParallelTest() {

	float refArr[8] = {1, 2, 3, 4, 5, 4, 3, 2};
	float nullArr[8] = {0, 0, 0, 0, 0, 0, 0, 0};

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[8];

	AVXUtils::cosineSimilarityIntraAVXParallel(
		refArr,
		nullArr,
		refArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		8,
		resultArr
		);

	for (int i = 0; i < 8; ++i) {
		if (i == 1) {
			QCOMPARE(MathUtils::tSame(resultArr[i], 1.0f), true);
			continue;
		}

		QCOMPARE(MathUtils::tZero(resultArr[i]), true);
	}

	float c[8] = {1, 1, 1, 1, 0, 0, 0, 0};
	float d[8] = {2, 0, 0, 0, 1, 1, 0, 0};
	AVXUtils::cosineSimilarityIntraAVXParallel(
		c,
		d,
		refArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		nullArr,
		8,
		resultArr
		);

	QVERIFY(MathUtils::tSame(resultArr[0], 0.408248f));
}

void AVXUtilsTests::sumParallelTest() {

	ERR_INIT

	float a[4] = {1, 1, 1, 1};
	float b[4] = {2, 2, 2, 2};
	float z[4] = {0, 0, 0, 0};

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[8];
	e = AVXUtils::sumParallel(
		4,
		a, b, z, z,z,z,z,z, resultArr
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(resultArr[0]), 4);
	QCOMPARE(static_cast<int>(resultArr[1]), 8);

}

void AVXUtilsTests::meanParallelTest() {

	ERR_INIT

	float a[4] = {1, 0, 1, 1};
	float b[4] = {2, 2, 2, 2};
	float z[4] = {0, 0, 0, 0};

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[8];
	e = AVXUtils::meanParallel(
		4,
		true,
		a, b, z, z,z,z,z,z, resultArr
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(resultArr[0]), 1);
	QCOMPARE(static_cast<int>(resultArr[1]), 2);

	std::memset(resultArr, 0, sizeof(resultArr));
	e = AVXUtils::meanParallel(
		4,
		false,
		a, b, z, z,z,z,z,z, resultArr
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(resultArr[0] * 100), static_cast<int>(0.75 * 100));
	QCOMPARE(static_cast<int>(resultArr[1]), 2);
}

void AVXUtilsTests::stDevParallelTest() {

	ERR_INIT

	float a[4] = {1, 0, 1, 1};
	float b[4] = {2, 2, 2, 2};
	float z[4] = {0, 0, 0, 0};
	float m[4] = {1, 2, 1, 2};

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float means[8];
	e = AVXUtils::meanParallel(
		4,
		true,
		a, b, z, m,z,z,z,z, means
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(means[0]), 1);
	QCOMPARE(static_cast<int>(means[1]), 2);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float stDevs[8];
	e = AVXUtils::stDevParallel(
		4,
		true,
		means,
		a, b, z, m,z,z,z,z,
		stDevs
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(stDevs[0]), 0);
	QCOMPARE(static_cast<int>(stDevs[1]), 0);
	QCOMPARE(std::isnan(stDevs[2]), true);
	QCOMPARE(static_cast<int>(stDevs[3] * 10), 5);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float means2[8];
	e = AVXUtils::meanParallel(
		4,
		false,
		a, b, z, m,z,z,z,z, means2
		);
	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(means2[0] * 100), static_cast<int>(0.75 * 100));
	QCOMPARE(static_cast<int>(means2[1]), 2);

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float stDevs2[8];
	e = AVXUtils::stDevParallel(
		4,
		false,
		means2,
		a, b, z, m,z,z,z,z,
		stDevs2
		);

	QCOMPARE(e, eNoError);
	QCOMPARE(static_cast<int>(stDevs2[0] * 1e6), 433012);
	QCOMPARE(stDevs2[1], 0);
	QCOMPARE(stDevs2[2], 0);
	QCOMPARE(static_cast<int>(stDevs2[3] * 10), 5);

}

QTEST_MAIN(AVXUtilsTests)

#include "AVXUtilsTests.moc"
