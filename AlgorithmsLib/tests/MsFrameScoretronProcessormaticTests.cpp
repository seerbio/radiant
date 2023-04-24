#include "MsFrameScoretronProcessormatic.h"

#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
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

private:

    static PythiaParameters pythiaParameters() {


        PythiaParameters pythiaParameters;

        pythiaParameters.returnPSMTopN = 1;
        pythiaParameters.maxTandemPointCount = 2;
        pythiaParameters.ms2ExtractionWidthPPM = 12.0;
        pythiaParameters.precursorExtractionWindowThomsons = 1.0;
        pythiaParameters.chargeStateMin = 2;
        pythiaParameters.chargeStateMax = 3;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
        );

        return pythiaParameters;
    }

    QMap<PeptideStringWithMods, QVector<MS2Ion>> buildPredictions() {

//        const QVector<MS2Ion> peptideQXNDYVAK = {
//                MS2Ion(129.066,0.14,"b1"), MS2Ion(130.086,0.36,"y1-NH3"), MS2Ion(147.113,0.69,"y1"), MS2Ion(200.139,0.14,"y2-H2O"),
//                MS2Ion(201.123,0.04,"y2-NH3"), MS2Ion(214.155,0.62,"a2"), MS2Ion(218.15,1,"y2"), MS2Ion(224.139,0.13,"b2-H2O"),
//                MS2Ion(225.123,0.25,"b2-NH3"), MS2Ion(242.15,0.12,"b2"), MS2Ion(299.208,0.03,"y3-H2O"), MS2Ion(317.218,0.32,"y3"),
//                MS2Ion(338.182,0.04,"b3-H2O"), MS2Ion(339.166,0.12,"b3-NH3"), MS2Ion(356.193,0.04,"b3"), MS2Ion(453.209,0.42,"b4-H2O"),
//                MS2Ion(454.193,0.13,"b4-NH3"), MS2Ion(462.271,0.04,"y4-H2O"), MS2Ion(471.22,0.08,"b4"), MS2Ion(480.282,0.39,"y4"),
//                MS2Ion(577.298,0.03,"y5-H2O"), MS2Ion(595.309,0.25,"y5"), MS2Ion(616.273,0.45,"b5-H2O"), MS2Ion(617.257,0.06,"b5-NH3"),
//                MS2Ion(634.283,0.03,"b5"), MS2Ion(691.341,0.08,"y6-H2O"), MS2Ion(692.325,0.21,"y6-NH3"), MS2Ion(709.352,0.61,"y6"),
//                MS2Ion(715.341,0.26,"b6-H2O"), MS2Ion(786.378,0.11,"b7-H2O"), MS2Ion(822.436,0.07,"y7"), MS2Ion(932.484,0.04,"b8-H2O"),
//                MS2Ion(950.494,0.02,"y8")
//        };
//
//        const QVector<MS2Ion> peptideSXSAAXEHK = {
//                MS2Ion(147.113, 0.05, "y1"), MS2Ion(173.128, 1, "a2"), MS2Ion(201.123, 0.2, "b2"),
//                MS2Ion(266.161, 0.08, "y2-H2O"),
//                MS2Ion(270.145, 0.04, "b3-H2O"), MS2Ion(284.172, 0.15, "y2"), MS2Ion(288.155, 0.05, "b3"),
//                MS2Ion(341.182, 0.02, "b4-H2O"),
//                MS2Ion(378.206, 0.02, "y7^2"), MS2Ion(395.204, 0.03, "y3-H2O"), MS2Ion(412.219, 0.02, "b5-H2O"),
//                MS2Ion(413.214, 0.08, "y3"),
//                MS2Ion(526.298, 0.06, "y4"), MS2Ion(597.335, 0.13, "y5"), MS2Ion(668.373, 0.07, "y6"),
//                MS2Ion(737.394, 0.02, "y7-H2O"),
//                MS2Ion(755.405, 0.39, "y7")
//        };
//
        QMap<PeptideStringWithMods, QVector<MS2Ion>> framePredictions = {
//                {"QXNDYVAK", peptideQXNDYVAK},
//                {"SXSAAXEHK", peptideSXSAAXEHK}
        };

        return framePredictions;
    }

};


void MsFrameScoretronProcessormaticTests::rescoreMsFrameTest() {

    ERR_INIT

    const QString scoreVectorsFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScores");

    const QString scansFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.prq.474966.frameScans");

    const int topN = 10;

    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> topCansInFrameIndex;
    MsFrameScoretronProcessormatic::processLogicForFrameScores(
            scoreVectorsFilePath,
            scansFilePath,
            topN,
            &topCansInFrameIndex
            );

}


QTEST_MAIN(MsFrameScoretronProcessormaticTests)

#include "MsFrameScoretronProcessormaticTests.moc"
