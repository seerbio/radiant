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

};


void AVXUtilsTests::copyAVXTest() {

	ERR_INIT

	size_t size = 96;
	alignas(32) float src[size] = {0};
	alignas(32) float dst[size] = {0};

	for (size_t i = 0; i < size; i++) {
		src[i] = static_cast<float>(i);
	}

	e = AVXUtils::copyAVXFFloat(src, dst, size);
	QCOMPARE(e, Err::eNoError);





}


QTEST_MAIN(AVXUtilsTests)

#include "AVXUtilsTests.moc"
