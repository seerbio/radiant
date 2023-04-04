//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"

#include "ParquetReaderBase.h"

#include <QtTest/QtTest>

class ParquetReaderBaseTests : public QObject
{
    Q_OBJECT

public:
    ParquetReaderBaseTests() = default;
    ~ParquetReaderBaseTests() override = default;

private Q_SLOTS:

    void qVectorToStdStringAndbytesStdStringToQVectorTest();
};

void ParquetReaderBaseTests::qVectorToStdStringAndbytesStdStringToQVectorTest() {

    const QVector<double> vecs = {666.6, 333.3, 666.6};

    const std::string bytesStdString = ParquetReaderBase::qVectorToBytesStdString(vecs);

    const QVector<double> vecDecoded = ParquetReaderBase::bytesStdStringToQVector<double>(bytesStdString);

    QCOMPARE(vecDecoded.size(), vecs.size());
    QCOMPARE(vecDecoded, vecs);
}

QTEST_MAIN(ParquetReaderBaseTests)
#include "ParquetReaderBaseTests.moc"
