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

    const QStringList cTermExpected = {"K", "R"};
    const QStringList nTermExpected = {};

    QCOMPARE(pythiaParameters.ms2ExtractionWidthPPM, 20.5);
    QCOMPARE(pythiaParameters.mzMinMS2, 100.1);
    QCOMPARE(pythiaParameters.mzMaxMS2, 1500.2);

    QCOMPARE(pythiaParameters.chargeStateMin, 1);
    QCOMPARE(pythiaParameters.chargeStateMax, 4);
    QCOMPARE(pythiaParameters.precursorExtractionWindowThomsons, 0.5);

    QCOMPARE(pythiaParameters.filterLength, 3);
    QCOMPARE(pythiaParameters.sigma, 1.0);
    QCOMPARE(pythiaParameters.signalToNoiseRatio, 2.0);
    QCOMPARE(pythiaParameters.smoothCount, 1);

    QCOMPARE(pythiaParameters.percentFDR, 2.0);
    QCOMPARE(pythiaParameters.reportDecoys, true);

}


QTEST_MAIN(PythiaParameterReaderTests)
#include "PythiaParameterReaderTests.moc"
