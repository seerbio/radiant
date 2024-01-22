#include "IsotopicDistributionBuilder.h"

#include "AminoAcids.h"
#include "MolecularFormula.h"

#include <QDebug>

#include <cmath>

namespace{

    QVector<double> normalizeToMax(const QVector<double> vec)
    {
        const double max = *std::max_element(vec.begin(), vec.end());

        QVector<double> qNormalized;
        std::transform(vec.begin(), vec.end(), std::back_inserter(qNormalized),
                       [max](double i) {return i / max;});

        return qNormalized;
    }


    double calculateMonoisotopeProbability(const MolecularFormula &molecularFormulaOfSequence) {

        return std::pow(C12_ABUNDANCE, molecularFormulaOfSequence.carbonCount)
               * std::pow(H1_ABUNDANCE , molecularFormulaOfSequence.hydrogenCount)
               * std::pow(N14_ABUNDANCE , molecularFormulaOfSequence.nitrogenCount)
               * std::pow(O16_ABUNDANCE , molecularFormulaOfSequence.oxygenCount)
               * std::pow(S32_ABUNDANCE , molecularFormulaOfSequence.sulfurCount);
    }


    QVector<double> calculateIsotopicDistribution(
            const MolecularFormula &molecularFormulaOfSequence,
            double stopThreshold = 0.001
                    ) {
        const double rootCarbon = -C12_ABUNDANCE / C13_ABUNDANCE;
        const double rootHydrogen = -H1_ABUNDANCE / H2_ABUNDANCE;
        const double rootNitrogen = -N14_ABUNDANCE / N15_ABUNDANCE;

        const double r1 = -O17_ABUNDANCE / (2 * O18_ABUNDANCE);
        const double r2 = std::sqrt(std::abs(std::pow(O17_ABUNDANCE,2) - (4*O16_ABUNDANCE*O18_ABUNDANCE))) / (2*O18_ABUNDANCE);
        const double rootOxygenModulus = std::sqrt(std::pow(r1, 2) + std::pow(r2, 2));
        const double rootOxygenArgument = std::atan(-r2/r1);

        const double firstPart = 0.000011149595585434016;
        const double s11 = -firstPart;
        const double s1i = -4.015974988771853;
        const double rootSulfurModulus1 = std::sqrt(std::pow(s11, 2) + std::pow(s1i, 2));
        const double rootSulferArgument1 = std::atan(s1i/s11);

        const double s21 = firstPart;
        const double s2i = -7.93332391982147;
        const double rootSulfurModulus2 = std::sqrt(std::pow(s21, 2) + std::pow(s2i, 2));
        const double rootSulferArgument2 = std::atan(s2i/s21);

        const double monoisotopeProb = calculateMonoisotopeProbability(molecularFormulaOfSequence);

        QVector<double> f;
        QVector<double> q;
        q.push_back(monoisotopeProb);

        int i = 1;
        while(q.back() > stopThreshold){

            const double fi = molecularFormulaOfSequence.carbonCount * std::pow(rootCarbon, -i)
                    + molecularFormulaOfSequence.hydrogenCount * std::pow(rootHydrogen, -i)
                    + molecularFormulaOfSequence.nitrogenCount * std::pow(rootNitrogen, -i)
                    + molecularFormulaOfSequence.oxygenCount * std::pow(std::abs(rootOxygenModulus), -i) * std::cos(-i * rootOxygenArgument)
                    + molecularFormulaOfSequence.sulfurCount * std::pow(std::abs(rootSulfurModulus1), -i) * std::cos(-i * rootSulferArgument1)
                    + molecularFormulaOfSequence.sulfurCount * std::pow(std::abs(rootSulfurModulus2), -i) * std::cos(-i * rootSulferArgument2)
                    ;

            double qI =  q[0] * fi;
            for(int j = 1; j < q.size(); j++){
                qI += q[j] * f[j-1];
            }

            f.push_front(fi);
            q.push_back(-qI/i);
            ++i;
        }

        return normalizeToMax(q);
    }

}//NAMESPACE
QVector<double> IsotopicDistributionBuilder::getIsotopicDistribution(const QStringList &sequence) {

    MolecularFormula molecularFormulaOfSequence
            = getMolecularFormulaOfBioPolymer(sequence);

    QVector<double> isotopicDistribution
            = calculateIsotopicDistribution(molecularFormulaOfSequence);

    return isotopicDistribution;
}


MolecularFormula
IsotopicDistributionBuilder::getMolecularFormulaOfBioPolymer(const QStringList &sequence) {

    const QMap<QChar, Molecule> aminoAcids = AminoAcids::aminoAcids();

    MolecularFormula returnFormula;
    for(const QChar &residue : sequence.join("")){

        Molecule currentFormula = aminoAcids.value(residue);
        returnFormula = returnFormula + currentFormula.molecularFormula();

    }

    return returnFormula + waterFormula;
}



namespace  {

    MolecularFormula convertPeptideMassToMolecularFormula(double mass) {

        AVERAGINE avg;
        const double numberOfAveragineUnits = mass / avg.averagineMass;

        MolecularFormula molForm;
        molForm.carbonCount = static_cast<int>(std::round(avg.C * numberOfAveragineUnits));
        molForm.hydrogenCount = static_cast<int>(std::round(avg.H * numberOfAveragineUnits));
        molForm.nitrogenCount = static_cast<int>(std::round(avg.N * numberOfAveragineUnits));
        molForm.oxygenCount = static_cast<int>(std::round(avg.O * numberOfAveragineUnits));
        molForm.sulfurCount = static_cast<int>(std::round(avg.S * numberOfAveragineUnits));

        const double averageMass = molForm.carbonCount * AverageMass::CARBON
                        + molForm.hydrogenCount * AverageMass::HYDROGEN
                        + molForm.nitrogenCount * AverageMass::NITROGEN
                        + molForm.oxygenCount * AverageMass::OXYGEN
                        + molForm.sulfurCount * AverageMass::SULFUR;

        const int hydrogenAdjustment = static_cast<int>(std::round(mass - averageMass));
        molForm.hydrogenCount += hydrogenAdjustment;

        return molForm;
    }

}//NAMESPACE
QVector<double> IsotopicDistributionBuilder::getIsotopicDistribution(double mass)
{

    MolecularFormula composition
            = convertPeptideMassToMolecularFormula(mass);

    return calculateIsotopicDistribution(composition);;
}




