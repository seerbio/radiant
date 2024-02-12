//
// Created by anichols on 2/5/24.
//

#ifndef PYTHIADIACPP_DEMULTIPLEXSCANERATOR_H
#define PYTHIADIACPP_DEMULTIPLEXSCANERATOR_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class MsScanInfo;

class ALGORITHMSFFLIB_EXPORTS DeMultiplexScanerator {

public:

    friend class DeMultiplexScaneratorTests;

    DeMultiplexScanerator(
            double ppmTol,
            double intensityFractionThreshold
            );

    ~DeMultiplexScanerator();

    Err deMultiplexScans(
            const QVector<ScanPoints*> &scans,
            const QVector<MsScanInfo> &msScanInfos,
            QVector<ScanPoints> *demultiplexedScans,
            QVector<MzTargetKey> *mzTargetKeys,
            double *windowSize
            );

private:

    Err _buildScanMaskMatrixTestAccess(
            const QVector<MsScanInfo> &msScanInfos,
            QVector<QVector<float>> *scanMaskMatrixVecs
            );

    Err _buildTransitionMatrixTestAccess(
            const QVector<ScanPoints*> &scans,
            QVector<QVector<float>> *transitionMatrixVecs
            );

private:

    Q_DISABLE_COPY(DeMultiplexScanerator) class Private;
    const QScopedPointer<Private> d_ptr;



};


#endif //PYTHIADIACPP_DEMULTIPLEXSCANERATOR_H
