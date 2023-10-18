#include "TargetDecoyCandidatePairManager.h"

#include <QtTest/QtTest>


class TargetDecoyCandidatePairManagerTests : public QObject
{
    Q_OBJECT

public:
    TargetDecoyCandidatePairManagerTests() = default;
    ~TargetDecoyCandidatePairManagerTests() override = default;

private Q_SLOTS:
    void loadModelTest();


};

void TargetDecoyCandidatePairManagerTests::loadModelTest() {

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


}



QTEST_MAIN(TargetDecoyCandidatePairManagerTests)
#include "TargetDecoyCandidatePairManagerTests.moc"
