//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "ErrorUtils.h"
#include "MsFraggerTronResultsReader.h"
#include "MsFraggerTronWorkFlow.h"
#include "MsReaderBase.h"
#include "ParallelUtils.h"


Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &pepLibUri
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(fragLibUri); ree;
    e = ErrorUtils::isNotEmpty(pepLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;
    m_pepLibUri = pepLibUri;

    ERR_RETURN

}

Err PythiaDIAWorkflow::processFile(const QString &mzmlFilePath) {

    ERR_INIT

    QString firstPassPSMsFilePath;
    QVector<TandemScanIon> tandemScanIons;
    e = runFirstPassMsFraggerTronWorkFlow(
            mzmlFilePath,
            &firstPassPSMsFilePath,
            &tandemScanIons
            ); ree;

    e = initReCalibratomatic(firstPassPSMsFilePath); ree;
    e = recalibrateTandemScanIons(&tandemScanIons); ree;
    e = optimizePythiaParameters(); ree;

    const QString outputFilePath
        = mzmlFilePath + ".recal" + S_GLOBAL_SETTINGS.DOT_PSM + S_GLOBAL_SETTINGS.DOT_CSV;

    e = runSecondPassMsFraggerTronWorkFlow(
            tandemScanIons,
            outputFilePath
            ); ree;

    qDebug() << "final output psm filepath" << outputFilePath;
    ERR_RETURN
}

Err PythiaDIAWorkflow::runFirstPassMsFraggerTronWorkFlow(
        const QString &mzmlFilePath,
        QString *firstPassPSMsFilePath,
        QVector<TandemScanIon> *tandemScanIons
        ) {

    ERR_INIT

    MsFraggerTronWorkFlow msFraggerTronWorkFlow;
    e = msFraggerTronWorkFlow.init(
            m_pythiaParameters,
            m_fragLibUri,
            m_pepLibUri
    );ree;

    e = msFraggerTronWorkFlow.processFile(
            mzmlFilePath,
            firstPassPSMsFilePath,
            tandemScanIons
    ); ree;

    ERR_RETURN
}

namespace {

    Err buildRecalibrationData(
            const QString &firstPassPSMCsvFilePath,
            QVector<InputSVM> *data
    ) {

        ERR_INIT

        QVector<RowToWrite> rowsToWrite;
        e = MsFraggerTronResultsReader::readCsv(
                firstPassPSMCsvFilePath,
                &rowsToWrite
        ); ree;



        int counter = 0;
        for (const RowToWrite &rtw : rowsToWrite) {

            counter++;

            if (rtw.isDecoy) {
                break;
            } //TODO think of a better stopping point.

            const QVector<double> &scanIonMZs = rtw.scanIonMZs;
            const QVector<double> &theoFragIonMZs = rtw.theoFragIonMZs;

            e = ErrorUtils::isEqual(scanIonMZs.size(), theoFragIonMZs.size()); ree;

            for (int i = 0; i < scanIonMZs.size(); i++) {

                const double scanMz = scanIonMZs.at(i);
                const double theoMz = theoFragIonMZs.at(i);
                const double ppmDiff = 1e6 * (scanMz - theoMz) / theoMz;

                InputSVM is;
                is.scanNumber = rtw.scanNumber;
                is.mzScan = scanMz;
                is.mzTheo = theoMz;
                is.ppmDiff = ppmDiff;

                data->push_back(is);
            }
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::initReCalibratomatic(const QString &firstPassPSMsFilePath) {

    ERR_INIT

    QVector<InputSVM> recalibrationData;
    e = buildRecalibrationData(
            firstPassPSMsFilePath,
            &recalibrationData
    ); ree;

    e = m_reCalibratomatic.initSVM(recalibrationData); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<TandemScanIon>> reCalParallelLogic(
            const QVector<TandemScanIon> &tandemScanIons,
            ReCalibratomatic reCalibratomatic
            ) {

        ERR_INIT

        QVector<TandemScanIon> reCalTandemScanIons;
        for (const TandemScanIon &tsi : tandemScanIons) {

            TandemScanIon reCalTsi = tsi;


            e = reCalibratomatic.recalibrateMz(
                    tsi.scanNumber,
                    tsi.mz,
                    &reCalTsi.mz
            );
            if (e != eNoError) {
                return {e, {}};
            }

            reCalTandemScanIons.push_back(reCalTsi);
        }

        return {e, reCalTandemScanIons};
    }

}
Err PythiaDIAWorkflow::recalibrateTandemScanIons(QVector<TandemScanIon> *tandemScanIons) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*tandemScanIons); ree;
    e = ErrorUtils::isTrue(m_reCalibratomatic.isInit()); ree;

    const int trancheBuffer = 0;

    QVector<QVector<TandemScanIon>> tranchedTandemScanIons;
    e = ParallelUtils::tranchVectorForParallelizationInOrder(
            *tandemScanIons
            , m_pythiaParameters.trancheSize,
            trancheBuffer,
            &tranchedTandemScanIons
            ); ree;

    const auto reCalBinder = std::bind(
            reCalParallelLogic,
            std::placeholders::_1,
            m_reCalibratomatic
    );

    QFuture<QPair<Err, QVector<TandemScanIon>>> futures = QtConcurrent::mapped(
            tranchedTandemScanIons,
            reCalBinder
            );
    futures.waitForFinished();

    tandemScanIons->clear();
    for (const QPair<Err, QVector<TandemScanIon>> &future : futures) {

        e = future.first; ree;
        tandemScanIons->append(future.second);
    }

    ERR_RETURN
}

Err PythiaDIAWorkflow::optimizePythiaParameters() {

    ERR_INIT

    const QPair<double, double> oldVsNewPPMStDev = m_reCalibratomatic.oldVsNewPPMStDev();
    const double oldPPMStDev = oldVsNewPPMStDev.first;
    const double newPPMStDev = oldVsNewPPMStDev.second;

    e = ErrorUtils::isTrue(oldPPMStDev > newPPMStDev); ree;
    e = ErrorUtils::isNotEqual(oldPPMStDev, 0.0); ree;

    const double ppmAdjustmentPostCal = newPPMStDev / oldPPMStDev;

    const double oldMs2ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPM;
    m_pythiaParameters.ms2ExtractionWidthPPM *= ppmAdjustmentPostCal;
    m_pythiaParameters.returnPSMTopN = 20; //TODO correct this.

    qDebug() << "PPM Extraction value adjusted from:" << oldMs2ExtractionWidthPPM
                << "to" << m_pythiaParameters.ms2ExtractionWidthPPM;

    ERR_RETURN
}

Err PythiaDIAWorkflow::runSecondPassMsFraggerTronWorkFlow(
        const QVector<TandemScanIon> &tandemScanIons,
        const QString &psmOutputFilePath
        ) {

    ERR_INIT

    MsFraggerTronWorkFlow msFraggerTronWorkFlow;
    e = msFraggerTronWorkFlow.init(
            m_pythiaParameters,
            m_fragLibUri,
            m_pepLibUri
    );ree;

    e = msFraggerTronWorkFlow.processScanIons(
            tandemScanIons,
            psmOutputFilePath
            ); ree;

    ERR_RETURN
}

