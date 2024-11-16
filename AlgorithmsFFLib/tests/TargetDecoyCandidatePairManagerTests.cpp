#include "TargetDecoyCandidatePairManager.h"

#include "FragLibReader.h"
#include "PythiaParameterReader.h"

#include <QtTest/QtTest>


class TargetDecoyCandidatePairManagerTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairManagerTests() = default;
    ~TargetDecoyCandidatePairManagerTests() override = default;

private Q_SLOTS:

    void initTest();
    void peptideStringWithModsFromPeptideSequenceChargeKeyTest();

};

void TargetDecoyCandidatePairManagerTests::initTest() {

    ERR_INIT

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            testFilePath,
            &fragLibReaderRows
    );
    QCOMPARE(e, eNoError);

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            &fragLibReaderRows
            );
    QCOMPARE(e, eNoError);

    fragLibReaderRows.clear();
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            &fragLibReaderRows
    );
    QCOMPARE(e, eError);

}

void TargetDecoyCandidatePairManagerTests::peptideStringWithModsFromPeptideSequenceChargeKeyTest() {

    ERR_INIT

    PeptideStringWithMods peptideStringWithMods;
    int charge;

    e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
            "BELLAFLOPSANDCHAUNCY|2",
            &peptideStringWithMods,
            &charge
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideStringWithMods, "BELLAFLOPSANDCHAUNCY");
    QCOMPARE(charge, 2);

    e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
            "BELLAFLOPSANDCHAUNCY",
            &peptideStringWithMods,
            &charge
    );
    QCOMPARE(e, eError);

    e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
            "BELLAFLOPSANDCHAUNCY||2",
            &peptideStringWithMods,
            &charge
    );
    QCOMPARE(e, eNoError);

    e = TargetDecoyCandidatePairManager::peptideStringWithModsFromPeptideSequenceChargeKey(
            "BELLAFLOPSANDCHAUNCY|R",
            &peptideStringWithMods,
            &charge
    );
    QCOMPARE(e, eError);

}


QTEST_MAIN(TargetDecoyCandidatePairManagerTests)
#include "TargetDecoyCandidatePairManagerTests.moc"
