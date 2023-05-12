//
// Created by anichols on 5/12/23.
//

#ifndef PYTHIADIACPP_EXTRACTEDSCANSREADER_H
#define PYTHIADIACPP_EXTRACTEDSCANSREADER_H

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

#include <QMap>
#include <vector>

using namespace Error;

namespace ExtractedScansReaderRowNamespace {

    const QString PEPTIDE_W_MODS = QStringLiteral("peptideWithMods");
    const QString MZ_FOUND = QStringLiteral("mzFound");
    const QString MZ_SEARCHED = QStringLiteral("mzSearched");
    const QString INTENSITY_FOUND = QStringLiteral("intensityFound");
    const QString INTENSITY_SEARCHED = QStringLiteral("intensitySearched");

    const QStringList keysToCheck = {
            PEPTIDE_W_MODS,
            MZ_FOUND,
            MZ_SEARCHED,
            INTENSITY_FOUND,
            INTENSITY_SEARCHED
    };
}

struct FILEREADERSLIB_EXPORTS ExtractedScansReaderRow : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    QVector<double> mzFound;
    QVector<double> mzSearched;
    QVector<double> intensityFound;
    QVector<double> intensitySearched;

    QMap<QString, QVariant> map() override {

        using namespace ExtractedScansReaderRowNamespace;

        return {
                {PEPTIDE_W_MODS, QVariant(peptideStringWithMods)},
                {MZ_FOUND, QVariant(qVectorToQByteArray(mzFound))},
                {MZ_SEARCHED, QVariant(qVectorToQByteArray(mzSearched))},
                {INTENSITY_FOUND, QVariant(qVectorToQByteArray(intensityFound))},
                {INTENSITY_SEARCHED, QVariant(qVectorToQByteArray(intensitySearched))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace ExtractedScansReaderRowNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideStringWithMods = dataMap.value(PEPTIDE_W_MODS).toString();
        mzFound = bytesArrayToQVector<double>(dataMap.value(MZ_FOUND).toByteArray());
        mzSearched = bytesArrayToQVector<double>(dataMap.value(MZ_SEARCHED).toByteArray());
        intensityFound = bytesArrayToQVector<double>(dataMap.value(INTENSITY_FOUND).toByteArray());
        intensitySearched = bytesArrayToQVector<double>(dataMap.value(INTENSITY_SEARCHED).toByteArray());

        ERR_RETURN
    }

};


#endif //PYTHIADIACPP_EXTRACTEDSCANSREADER_H
