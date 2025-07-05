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


QTEST_MAIN(AVXUtilsTests)

#include "AVXUtilsTests.moc"
