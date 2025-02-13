//
// Created by andrewnichols on 2/11/25.
//

#ifndef OPTIMIZORMSCALIBRATOMATIC_H
#define OPTIMIZORMSCALIBRATOMATIC_H

#include "AlgorithmsOptimizersLib_Exports.h"
#include "GlobalSettings.h"


using namespace Error;


class ALGORITHMSOPTIMIZERSLIB_EXPORTS OptimizorMsCalibratomatic {

public:
    OptimizorMsCalibratomatic() = default;
    ~OptimizorMsCalibratomatic() = default;

    static Err optimize(
        const QVector<QString> &msDataFilePaths,
        const QString &libFilePath,
        const QString &parametersFilePath,
        int maxGenerationIterations
        );

private:

};



#endif //OPTIMIZORMSCALIBRATOMATIC_H
