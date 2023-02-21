//
// Created by anichols on 2/19/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYTRON_H
#define PYTHIADIACPP_FRAGLIBRARYTRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "PythiaParameterReader.h"

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

};

struct ALGORITHMSLIB_EXPORTS FragLibIonTranche {
    QVector<FragLibIon> fragLibIonsTranche;
    double mzMin = -1;
    double mzMax = -1;
};


class ALGORITHMSLIB_EXPORTS FragLibraryTron {

public:

    friend class FragLibraryTronTests;

    explicit FragLibraryTron(const PythiaParameters &pythiaParameters);
    ~FragLibraryTron() = default;



    static Err writeFragLibIons(
            const QVector<FragLibIon> &fragLibIons,
            const QString &fragLibIonsFilePath
            );

    Err readFragLibIons(const QString &fragLibIonsFilePath);

    Err getFragLibIonTranches(
            const QVector<QPair<MzMin, MzMax>> &tranchLimits,
            QVector<FragLibIonTranche> *fragLibIonTranches
            );

    QVector<FragLibIon> fragLibIons();

private:

    QPair<int, int> getFragLibIonTrancheStartStop(const QPair<MzMin, MzMax> &tranchLimits);

private:

    PythiaParameters m_pythiaParameters;
    QVector<FragLibIon> m_fragLibIons;

};


#endif //PYTHIADIACPP_FRAGLIBRARYTRON_H
