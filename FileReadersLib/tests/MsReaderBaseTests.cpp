//
// Created by anichols on 11/07/2021.
//

#include "MsReaderBase.h"
#include "MsReaderPointerFactory.h"

#include <QtTest/QtTest>


class MsReaderBaseTests : public QObject
{
    Q_OBJECT

public:
    MsReaderBaseTests() = default;
    ~MsReaderBaseTests() override = default;

private Q_SLOTS:

    void openMzMlTest();
    void openPrqTest();
    void isDIATest();


private:
    //TODO use proper path procedures.
    const QString m_mzMLFilePath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

    //TODO use proper path procedures.
    const QString m_prqFilePath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML.prq");

};

void MsReaderBaseTests::openMzMlTest() {

    QSKIP("SKIPPING FOR DEV");

    QPair<Err, MsReaderPointer> msReaderFactoryResult = MsReaderPointerFactory::createInstance(m_mzMLFilePath);
    QCOMPARE(msReaderFactoryResult.first, eNoError);

}

void MsReaderBaseTests::openPrqTest() {

    QPair<Err, MsReaderPointer> msReaderFactoryResult = MsReaderPointerFactory::createInstance(m_prqFilePath);
    QCOMPARE(msReaderFactoryResult.first, eNoError);
}

void MsReaderBaseTests::isDIATest() {

    QPair<Err, MsReaderPointer> msReaderFactoryResult = MsReaderPointerFactory::createInstance(m_prqFilePath);
    QCOMPARE(msReaderFactoryResult.first, eNoError);

    MsReaderPointer msReaderPointer = msReaderFactoryResult.second;

    const bool msFileIsDIA = msReaderPointer->isDIA();
    QCOMPARE(msFileIsDIA, true);
}


QTEST_MAIN(MsReaderBaseTests)
#include "MsReaderBaseTests.moc"
