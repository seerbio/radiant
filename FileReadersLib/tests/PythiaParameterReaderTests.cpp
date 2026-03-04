//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "PythiaParameterReader.h"

#include <QtTest/QtTest>
#include <QVariant>

class PythiaParameterReaderTests : public QObject
{
    Q_OBJECT

public:
    PythiaParameterReaderTests() = default;
    ~PythiaParameterReaderTests() override = default;

private Q_SLOTS:

    static void readFileTest();

};

void PythiaParameterReaderTests::readFileTest() {

    ERR_INIT

    const QString filePath
            = QDir(qApp->applicationDirPath()).filePath("test_params_wide_window.pythiaConfig");

    PythiaParameterReader reader;

    PythiaParameters pythiaParameters;
    PythiaParameterReader::buildPythiaParameters(filePath, &pythiaParameters);
    QCOMPARE(e, eNoError);

    QCOMPARE(pythiaParameters.threadCount, 16);
    QCOMPARE(pythiaParameters.verbosity, 1);
    QCOMPARE(pythiaParameters.writeRadiantDIA, true);
    QCOMPARE(pythiaParameters.chargeStateMin, 1);
    QCOMPARE(pythiaParameters.chargeStateMax, 4);
    QCOMPARE(pythiaParameters.mzMinMS2, 200.0);
    QCOMPARE(pythiaParameters.mzMaxMS2, 1500.0);
    QCOMPARE(pythiaParameters.peptideLengthMin, 8);
    QCOMPARE(pythiaParameters.peptideLengthMax, 31);
    QCOMPARE(pythiaParameters.trancheSizeMax, 1e4);
    QCOMPARE(pythiaParameters.precursorExtractionWindowThomsons, 0.5);
    QCOMPARE(pythiaParameters.ms1ExtractionWidthPPM, 20.0);
    QCOMPARE(pythiaParameters.filterLengthIntegration, 6);
    QCOMPARE(pythiaParameters.filterLengthMS2, 4);
    QCOMPARE(pythiaParameters.ionsSharedToReject, 2);
    QCOMPARE(pythiaParameters.ms2ExtractionWidthPPM, 21.0);
    QCOMPARE(pythiaParameters.minMs2FragCount, 3);
    QCOMPARE(pythiaParameters.scanTimeWindowStDevs, 4);
    QCOMPARE(pythiaParameters.subtractShadows, false);
    QCOMPARE(pythiaParameters.smoothCountMS2, 2);
    QCOMPARE(pythiaParameters.stopThresholdFractionMS2, 0.666f);
    QCOMPARE(pythiaParameters.percentFDR, 2.0);
    QCOMPARE(pythiaParameters.reportDecoys, true);
    QCOMPARE(pythiaParameters.filterLength, 4);
    QCOMPARE(pythiaParameters.sigma, 1.1);
    QCOMPARE(pythiaParameters.signalToNoiseRatio, 2.1);
    QCOMPARE(pythiaParameters.smoothCount, 3);
    QCOMPARE(pythiaParameters.minScanCount, 4);
    QCOMPARE(pythiaParameters.skipScanCount, 4);

    pythiaParameters.print();

}


QTEST_MAIN(PythiaParameterReaderTests)
#include "PythiaParameterReaderTests.moc"
