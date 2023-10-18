//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePairManager.h"

#include "FragLibReader.h"

namespace {

    Err buildMzIons(
            const FragLibReaderRow &row,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(row.mzVals); ree;
        e = ErrorUtils::isEqual(row.mzVals.size(), row.intensityVals.size()); ree;

        const QStringList ionLabels = row.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

        e = ErrorUtils::isEqual(row.mzVals.size(), ionLabels.size()); ree;

        for (int i = 0; i < row.mzVals.size(); i++) {
            MS2Ion ion;
            ion.mz = row.mzVals.at(i);
            ion.intensity = row.intensityVals.at(i);
            ion.ionLabel = ionLabels.at(i);
            ion.iRT = static_cast<float>(row.iRT);
            ion.charge = ion.ionLabel.contains("^2") ? 2 : 1;
            ms2Ions->push_back(ion);
        }

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairManager::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibFileUri
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibFileUri); ree;

    FragLibReader fragLibReader;
    e = fragLibReader.getMS2IonsTopN()


    ERR_RETURN
}
