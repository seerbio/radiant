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

    void isDIATest();

};


void MsReaderBaseTests::isDIATest() {


    ERR_INIT

    const QString &prqFFFilePath = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(prqFFFilePath);
    QCOMPARE(e, eNoError);

    const bool msFileIsDIA = msReaderPointerAcc.ptr->isDIA();
    QCOMPARE(msFileIsDIA, true);
}


QTEST_MAIN(MsReaderBaseTests)
#include "MsReaderBaseTests.moc"
