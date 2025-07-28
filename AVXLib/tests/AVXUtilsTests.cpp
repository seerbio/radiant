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

	const size_t size = 4;

	QVector<float> v1(size, 1);
	QVector<float> v2(size, 2);
	QVector<float> v3(size, 3);
	QVector<float> vAsc(size);
	std::iota(vAsc.begin(), vAsc.end(), 1.0f);

	const QVector<float> kernel = {1, 2, 1};

	e = AVXUtils::convolveWithKernelAVXFloat(
		kernel,
		size,
		vAsc.data(),
		v2.data(),
		v3.data(),
		v1.data(),
		v2.data(),
		v3.data(),
		v1.data(),
		v2.data()
		);

	const QVector<int> resultAsc = {4, 6, 9, 3};
	const QVector<int> resultV1 = {3, 4, 4, 3};

	for (int i = 0; i < size; i++) {
		QCOMPARE(static_cast<int>(vAsc[i]), resultAsc[i]);
		QCOMPARE(static_cast<int>(v1[i]), resultV1[i]);
	}

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

	AVXUtils::printAVXInt32(int32Output1);
	AVXUtils::printAVXInt32(int32Output2);

	__m256i int32OutputRecombine = _mm256_packs_epi32(int32Output1, int32Output2);
	AVXUtils::printAVXInt16(int32OutputRecombine);

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
	AVXUtils::subtractArraysAVX2(arr1, arr2, 16);

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

QTEST_MAIN(AVXUtilsTests)

#include "AVXUtilsTests.moc"
