//
// Created by anichols on 2/19/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYTRON_H
#define PYTHIADIACPP_FRAGLIBRARYTRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QDataStream>


using namespace Error;


struct ALGORITHMSLIB_EXPORTS FragLibIon {

        int peptideId = -1;
        double mzFrag = -1.0;
        double peptideMass = -1.0;

        FragLibIon() = default;

        FragLibIon(
                int peptideId,
                double mzFrag,
                double peptideMass
                )
                : peptideId(peptideId)
                , mzFrag(mzFrag)
                , peptideMass(peptideMass)
                {}

        friend QDataStream &operator <<(QDataStream &stream, const FragLibIon &fli) {
            stream << fli.peptideId;
            stream << fli.mzFrag;
            stream << fli.peptideMass;

            return stream;
        }

        friend QDataStream &operator >>(QDataStream &stream, FragLibIon &fli) {
            stream >> fli.peptideId;
            stream >> fli.mzFrag;
            stream >> fli.peptideMass;

            return stream;
        }

    [[nodiscard]] QString hashedKey() const {

        const int hashingPrecision = 4;

        const int mzHashed = MathUtils::hashDecimal(mzFrag, hashingPrecision);
        const int targetHashed = MathUtils::hashDecimal(peptideMass, hashingPrecision);

        return S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT
                .arg(mzHashed)
                .arg(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP)
                .arg(targetHashed);
    }

    static QPair<MZION , TARGETMZ> unhashKey(const QString &hashedKey) {
        const QStringList splitStr = hashedKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT,
                QString::SkipEmptyParts
        );

        const int hashingPrecision = 4;

        return {MathUtils::unHashDecimal<double>(splitStr.first().toInt(), hashingPrecision),
                MathUtils::unHashDecimal<double>(splitStr.back().toInt(), hashingPrecision)};
    }
};

class ALGORITHMSLIB_EXPORTS FragLibraryTron {

public:

    friend class FragLibraryTronTests;

    FragLibraryTron() = default;
    ~FragLibraryTron() = default;

    static Err writeFragLibIons(
            const QVector<FragLibIon> &fragLibIons,
            const QString &fragLibIonsFilePath
            );

    Err readFragLibIons(const QString &fragLibIonsFilePath);

    Err buildFragLibIonsUniqueIndexes();

private:

    QVector<FragLibIon> m_fragLibIons;
    QMap<MzTargetKey , QVector<FragLibIonPeptideId>> m_fragLibIonsUniquePeptideId;

};


#endif //PYTHIADIACPP_FRAGLIBRARYTRON_H
