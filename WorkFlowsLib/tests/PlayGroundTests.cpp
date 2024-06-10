//
// Created by anichols on 11/07/2021.
//

#include "../../UtilsLib/src/ErrorUtils.h"
#include "../../FileReadersLib/src/MsReaderBase.h"

#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "../../UtilsLib/src/GlobalSettings.h"
#include "../../FileReadersLib/src/MsReaderBase.h"
#include "../../FileReadersLib/src/MsReaderParquet.h"
#include "../../FileReadersLib/src/MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePairScoretron.h"


class PlayGroundTests : public QObject
{
    Q_OBJECT

public:
    PlayGroundTests() = default;
    ~PlayGroundTests() override = default;

private Q_SLOTS:

void testme();

};

void PlayGroundTests::testme() {

    ERR_INIT

    const ScanNumber scanNumber = 11114;
    const QString targetKey = "625034";

    const QString fragLibFilepath = QStringLiteral("/home/anichols/Desktop/Data/Libraries/predicted_lib.speclib.fragLibFF");
    const QString msDataFilepath = QStringLiteral("/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML.prqFF");

    MsReaderPointerAcc msReader;
    e = msReader.openFile(msDataFilepath);
    QCOMPARE(e, eNoError);



}

QTEST_MAIN(PlayGroundTests)
#include "PlayGroundTests.moc"
