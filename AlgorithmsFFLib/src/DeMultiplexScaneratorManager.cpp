//
// Created by anichols on 2/9/24.
//

#include "DeMultiplexScaneratorManager.h"

#include "DeMultiplexScanerator.h"
#include "MsReaderPointerAcc.h"

#include <QtConcurrent/QtConcurrent>

namespace {

    struct DemultiParallelInput {
        QVector<ScanPoints*> scans;
        QVector<MsScanInfo> msScanInfos;
    };

    struct DemultiParallelOutput {
        QVector<ScanPoints> demultiplexedScans;
        QVector<MzTargetKey> mzTargetKeys;
        double windwoSize = -1.0;
    };

    QPair<Err, DemultiParallelOutput> demultiplexLogic(const DemultiParallelInput &demultiParallelInput) {

        ERR_INIT

        DemultiParallelOutput demultiParallelOutput;
        DeMultiplexScanerator deMultiplexScanerator(10.0, 0.9);
        e = deMultiplexScanerator.deMultiplexScans(
                demultiParallelInput.scans,
                demultiParallelInput.msScanInfos,
                &demultiParallelOutput.demultiplexedScans,
                &demultiParallelOutput.mzTargetKeys,
                &demultiParallelOutput.windwoSize
        );

        return {e, demultiParallelOutput};
    }

}//namespace
Err DeMultiplexScaneratorManager::demultiplexMsReader(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    const int msLevel = 2;
    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc->ptr->getMsScanInfos(msLevel);

    const QHash<MzTargetKey, bool> lastTargetKeys = {
            {"506000", true}
    };

    QVector<QVector<MsScanInfo>> groupedMsScanInfosAll;

    QVector<MsScanInfo> groupedMsScanInfosCurrent;
    for(const MsScanInfo &msi : msScanInfos) {
        if (!lastTargetKeys.contains(msi.targetKey())) {
            groupedMsScanInfosCurrent.push_back(msi);
            continue;
        }

        groupedMsScanInfosCurrent.push_back(msi);
        groupedMsScanInfosAll.push_back(groupedMsScanInfosCurrent);
        groupedMsScanInfosCurrent.clear();
    }

    QMap<ScanNumber, ScanPoints*> scanPoints;
    e = msReaderPointerAcc->ptr->getScanPoints(2, &scanPoints); ree;

    QVector<QVector<ScanPoints*>> groupedScanPointsAll;
    QVector<ScanPoints*> groupedScanPointsCurrent;
    for (const QVector<MsScanInfo> &msis : groupedMsScanInfosAll) {

        for (const MsScanInfo &m : msis) {
            ScanPoints *sp = scanPoints.value(m.scanNumber);
            groupedScanPointsCurrent.push_back(sp);

        }

        groupedScanPointsAll.push_back(groupedScanPointsCurrent);
        groupedScanPointsCurrent.clear();
    }

    e = ErrorUtils::isEqual(groupedMsScanInfosAll.size(), groupedScanPointsAll.size()); ree;

    QVector<DemultiParallelInput> inputs;
    for (int i = 0; i < groupedScanPointsAll.size(); i++) {
        DemultiParallelInput demultiParallelInput;
        demultiParallelInput.msScanInfos = groupedMsScanInfosAll.at(i);
        demultiParallelInput.scans = groupedScanPointsAll.at(i);

        inputs.push_back(demultiParallelInput);
    }

    QFuture<QPair<Err, DemultiParallelOutput>> results = QtConcurrent::mapped(
            inputs,
            demultiplexLogic
    );
    results.waitForFinished();

    QMap<ScanNumber, ScanPoints> newScanPoints;
    QMap<ScanNumber, MsScanInfo> newMsScanInfos;

    int newScanNumber = 0;
    int outputCounter = 0;
    for (const QPair<Err, DemultiParallelOutput> &res : results) {

        e = res.first; ree;

        const DemultiParallelInput &demultiParallelInput = inputs.at(outputCounter++);
        const DemultiParallelOutput &demultiParallelOutput = res.second;

        e = ErrorUtils::isEqual(
                demultiParallelOutput.mzTargetKeys.size(),
                demultiParallelOutput.demultiplexedScans.size()
        ); ree;

        for (int i = 0; i < demultiParallelOutput.mzTargetKeys.size(); i++) {

            MsScanInfo newMsScanInfo;
            newMsScanInfo.isoWindowLower = demultiParallelOutput.windwoSize / 2.0;
            newMsScanInfo.isoWindowUpper = demultiParallelOutput.windwoSize / 2.0;
            newMsScanInfo.scanNumber = newScanNumber++;
            newMsScanInfo.msLevel = 2;
            newMsScanInfo.ionMobilityIndex = -1;
            newMsScanInfo.ionMobilityDriftTime = -1.0;
            newMsScanInfo.scanTime = demultiParallelInput.msScanInfos.front().scanTime;
            newMsScanInfo.collisionEnergy = demultiParallelInput.msScanInfos.front().collisionEnergy;
            e = ErrorUtils::toFloat(demultiParallelOutput.mzTargetKeys.at(i) , &newMsScanInfo.precursorTargetMz); ree;

            newMsScanInfos.insert(newMsScanInfo.scanNumber, newMsScanInfo);
            newScanPoints.insert(newMsScanInfo.scanNumber, demultiParallelOutput.demultiplexedScans.at(i));
        }

    }

    msReaderPointerAcc->ptr->setMsScanInfo(newMsScanInfos);
    e = msReaderPointerAcc->ptr->setScanPoints(newScanPoints); ree;

    ERR_RETURN
}
