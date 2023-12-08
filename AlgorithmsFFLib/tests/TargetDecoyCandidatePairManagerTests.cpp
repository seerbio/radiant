#include "TargetDecoyCandidatePairManager.h"

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
    void getTargetDecoyCandidatePairPointersTest1();
    void getTargetDecoyCandidatePairPointersTest2();
    void peptideStringWithModsFromPeptideSequenceChargeKeyTest();

};

void TargetDecoyCandidatePairManagerTests::initTest() {

    ERR_INIT

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibDF");

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
            );
    QCOMPARE(e, eNoError);

    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath + QStringLiteral("BelltrixAndKalliope")
    );
    QCOMPARE(e, eFileError);

}

void TargetDecoyCandidatePairManagerTests::getTargetDecoyCandidatePairPointersTest1() {

    ERR_INIT

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibDF");

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
    );
    QCOMPARE(e, eNoError);

    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    e = targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
            666.0,
            667.0,
            &targetDecoyPointers
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(targetDecoyPointers.size(), 25);

    std::sort(
            targetDecoyPointers.begin(),
            targetDecoyPointers.end(),
            [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){return l->mz() < r->mz();}
            );
    QVERIFY(MathUtils::tSame(targetDecoyPointers.front()->mz(), 666.016));
    QVERIFY(MathUtils::tSame(targetDecoyPointers.back()->mz(), 666.88));

}

void TargetDecoyCandidatePairManagerTests::getTargetDecoyCandidatePairPointersTest2() {

    ERR_INIT

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibDF");

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
    );
    QCOMPARE(e, eNoError);

    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    e = targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
            666.0,
            667.0,
            0.5,
            &targetDecoyPointers
    );
    QCOMPARE(e, eNoError);
    QVERIFY(targetDecoyPointers.size() < 25);

    std::sort(
            targetDecoyPointers.begin(),
            targetDecoyPointers.end(),
            [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){return l->mz() < r->mz();}
    );
    QVERIFY(targetDecoyPointers.front()->mz() >= 666.0);
    QVERIFY(targetDecoyPointers.back()->mz() <= 667.00);
}

void TargetDecoyCandidatePairManagerTests::peptideStringWithModsFromPeptideSequenceChargeKeyTest() {

    ERR_INIT

    PeptideStringWithMods peptideStringWithMods;
    Charge charge;

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
