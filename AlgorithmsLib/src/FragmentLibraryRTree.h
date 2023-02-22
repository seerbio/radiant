//
// Created by Drucifer on 7/9/2022.
//

#ifndef PYTHIACPP_FRAGMENTLIBRARYRTREE_H
#define PYTHIACPP_FRAGMENTLIBRARYRTREE_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibraryTron.h"
#include "GlobalSettings.h"

#include <QScopedPointer>


using namespace Error;


class ALGORITHMSLIB_EXPORTS FragmentLibraryRTree {


public:

    FragmentLibraryRTree();

    ~FragmentLibraryRTree();

    Err init(
            const QVector<FragLibIon> &fragLibIons,
            const QPair<int, int> &minMaxCharge,
            double ppmTolerance,
            double precursorExtractionWindowThomsons
            );

    QHash<PeptideId, MZION> getPeptidesTableIds(
            double mz,
            double targetMass,
            const QPair<double, double> &targetWindow
    );

    int size();

    void setKey(int key);

    int getKey();


private:

    Q_DISABLE_COPY(FragmentLibraryRTree) class Private;
    const QScopedPointer<Private> d_ptr;
    int m_key;

};


#endif //PYTHIACPP_FRAGMENTLIBRARYRTREE_H
