//
// Created by anichols on 7/1/23.
//

#ifndef PYTHIADIACPP_FRAGLIBIONRTREE_H
#define PYTHIADIACPP_FRAGLIBIONRTREE_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"


using namespace Error;

class MS2IonsSeparated;

struct FragLibIon {
    PeptideId peptideId = -1;
    double mz = -1.0;
    double iRT = -1.0;
    double iMobility = -1.0;
    double intensity = -1.0;
    IonIndex ionIndex = -1;
    IonType ionType;
    Charge charge = -1;
};


class ALGORITHMSLIB_EXPORTS FragLibIonRTree {

public:

    FragLibIonRTree();
    ~FragLibIonRTree();

    Err init(const QMap<PeptideStringWithMods, MS2IonsSeparated> &fragPreds);

private:

    Q_DISABLE_COPY(FragLibIonRTree) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_FRAGLIBIONRTREE_H
