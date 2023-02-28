//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "ErrorUtils.h"
#include "MsFraggerTronResultsReader.h"
#include "MsFraggerTronWorkFlow.h"


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
    e = runFirstPassMsFraggerTronWorkFlow(
            mzmlFilePath,
            &firstPassPSMsFilePath
            ); ree;

    e = initReCalibratomatic(firstPassPSMsFilePath); ree;


    ERR_RETURN
}

Err PythiaDIAWorkflow::runFirstPassMsFraggerTronWorkFlow(
        const QString &mzmlFilePath,
        QString *firstPassPSMsFilePath
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
            firstPassPSMsFilePath
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

