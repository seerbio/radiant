//
// Created by andrewnichols on 8/15/24.
//

#ifndef TURBOXIC2D_H
#define TURBOXIC2D_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "PointFF.h"
#include "TurboXIC.h"

class ALGORITHMSFFLIB_EXPORTS TurboXIC2D {

public:

    TurboXIC2D();
    ~TurboXIC2D();

    /**
    * @brief Initializes the TurboXIC private implementation.
    *
    * This method initializes the TurboXIC private implementation by populating the scan number vs. intensity pointers,
    * creating an R-tree for spatial indexing, and performing necessary checks.
    *
    * @param scanNumberVsScanPoints A map of scan numbers to scan points.
    * @return An error code indicating success or failure.
    */
    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) const;


    XICPoints extractPointsXIC(
            float dim1Min,
            float dim1Max,
            float dim2Min,
            float dim2Max
    );

    Err getRTreeLimits(
        float *dim1Min,
        float *dim1Max,
        float *dim2Min,
        float *dim2Max
        ) const;

    [[nodiscard]] bool isInit() const;

private:

    Q_DISABLE_COPY(TurboXIC2D) class Private;
    const QScopedPointer<Private> d_ptr;

};



#endif //TURBOXIC2D_H
