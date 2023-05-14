#ifndef ISOTOPICDISTRIBUTIONBUILDER_H
#define ISOTOPICDISTRIBUTIONBUILDER_H

#include "AlgorithmsLib_Exports.h"
#include "Molecule.h"

#include <QMap>
#include <QVector>
#include <QPoint>



class ALGORITHMSLIB_EXPORTS IsotopicDistributionBuilder {


public:

    static QVector<double> getIsotopicDistribution(const QStringList &sequence);

    static QVector<double> getIsotopicDistribution(double mass);

    static MolecularFormula getMolecularFormulaOfBioPolymer(const QStringList &sequence);

};

#endif // ISOTOPICDISTRIBUTIONBUILDER_H
