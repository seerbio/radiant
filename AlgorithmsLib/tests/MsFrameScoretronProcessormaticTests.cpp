#include "MsFrameScoretronProcessormatic.h"

#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsFrameScoretron.h"
#include "FrameExtractsReader.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"

#include <QElapsedTimer>
#include <QVector>
#include <QtTest>


class MsFrameScoretronProcessormaticTests : public QObject
{

    Q_OBJECT
    
public:
    MsFrameScoretronProcessormaticTests() = default;
    ~MsFrameScoretronProcessormaticTests() override = default;


private slots:

    void rescoreMsFrameTest();


};


void MsFrameScoretronProcessormaticTests::rescoreMsFrameTest() {

    ERR_INIT

    const QString frameExtractionFilepath
        = "/home/anichols/Desktop/Testing/EXP22092_2022ms0742X32_A.raw.mzML.reCal.prq.635039.frameExtractions";

    QVector<FrameExtractsReaderRow> frameExtracts;
    e = ParquetReader::read(
            frameExtractionFilepath,
            &frameExtracts
            );

    QCOMPARE(e, eNoError);

    for (const FrameExtractsReaderRow &fer : frameExtracts) {

        if (fer.frameIndexApex != 191) {
            continue;
        }

        if (fer.yIonMzStdActual.isEmpty()) {
            continue;
        }

        ScanPoints scanPoints;
        e = ParallelUtils::zip(
                fer.yIonMzValsActual,
                fer.yIonIntesitiesActual,
                &scanPoints
                );
        QCOMPARE(e, eNoError);

        qDebug() << scanPoints;
    }




}


QTEST_MAIN(MsFrameScoretronProcessormaticTests)

#include "MsFrameScoretronProcessormaticTests.moc"
