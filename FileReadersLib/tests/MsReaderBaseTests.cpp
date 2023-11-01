//
// Created by anichols on 11/07/2021.
//

#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"

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
    QSKIP("TODO: enable with internal test data");

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_mzMLFilePath);
    QCOMPARE(e, eNoError);

}

void MsReaderBaseTests::openPrqTest() {
    QSKIP("TODO: enable with internal test data");

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFilePath);
    QCOMPARE(e, eNoError);
}

void MsReaderBaseTests::isDIATest() {
    QSKIP("TODO: enable with internal test data");

    QSKIP("activate when proper pathing is used");

    ERR_INIT

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(m_prqFilePath);
    QCOMPARE(e, eNoError);

    const bool msFileIsDIA = msReaderPointerAcc.ptr->isDIA();
    QCOMPARE(msFileIsDIA, true);
}


QTEST_MAIN(MsReaderBaseTests)
#include "MsReaderBaseTests.moc"
