#include "TargetDecoyCandidatePairManager.h"

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


};

void TargetDecoyCandidatePairManagerTests::initTest() {

    ERR_INIT

    //TODO accumulate this test file to S3 for DL for unit test
    //TODO use smaller file to save on DL time.
    const QString testFilePath
            = "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta.predicted.speclib.fragLib";

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
            );
    QCOMPARE(e, eNoError);

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManagerParameterFail;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameters(),
            testFilePath
    );
    QCOMPARE(e, eError);

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManagerFilePathFail;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath + "KalliopeAndBellatrix"
    );
    QCOMPARE(e, eFileError);


}

void TargetDecoyCandidatePairManagerTests::getTargetDecoyCandidatePairPointersTest1() {

    ERR_INIT

    //TODO accumulate this test file to S3 for DL for unit test
    //TODO use smaller file to save on DL time.
    const QString testFilePath
            = "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta.predicted.speclib.fragLib";

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
    );
    QCOMPARE(e, eNoError);

    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    e = targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(700.0, 710.0, &targetDecoyPointers);
    QCOMPARE(e, eNoError);
    QCOMPARE(targetDecoyPointers.size(), 15792);

    TargetDecoyCandidatePair* tdcp = targetDecoyPointers.at(666);
    QVERIFY(MathUtils::tSame(tdcp->mz(), 700.381));
    QVERIFY(MathUtils::tSame(tdcp->mass(), 1398.75));
}

void TargetDecoyCandidatePairManagerTests::getTargetDecoyCandidatePairPointersTest2() {

    ERR_INIT

    //TODO accumulate this test file to S3 for DL for unit test
    //TODO use smaller file to save on DL time.
    const QString testFilePath
            = "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta.predicted.speclib.fragLib";

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFilePath
    );
    QCOMPARE(e, eNoError);

    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    e = targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(700.0, 710.0, 0.1, &targetDecoyPointers);
    QCOMPARE(e, eNoError);
    QCOMPARE(targetDecoyPointers.size(), 1579);

}


QTEST_MAIN(TargetDecoyCandidatePairManagerTests)
#include "TargetDecoyCandidatePairManagerTests.moc"
