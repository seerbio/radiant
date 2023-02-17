#ifndef MOLECULARFORMULAS_H
#define MOLECULARFORMULAS_H

#include "ChemLib_Exports.h"
#include "Error.h"

#include <QString>

using namespace Error;

namespace MolecularFormulas
{

    struct CHEMLIB_EXPORTS MolecularFormula{

        int carbonCount = 0;
        int hydrogenCount = 0;
        int nitrogenCount = 0;
        int oxygenCount = 0;
        int sulfurCount = 0;
        int phosphorusCount = 0;

        MolecularFormula() = default;

        MolecularFormula(int carbonCount,
                         int hydrogenCount,
                         int nitrogenCount,
                         int oxygenCount,
                         int sulfurCount,
                         int phosphorusCount)
                         : carbonCount(carbonCount)
                         , hydrogenCount(hydrogenCount)
                         , nitrogenCount(nitrogenCount)
                         , oxygenCount(oxygenCount)
                         , sulfurCount(sulfurCount)
                         , phosphorusCount(phosphorusCount){}

        MolecularFormula operator+(const MolecularFormula &a) {
            MolecularFormula c;
            c.carbonCount = this->carbonCount + a.carbonCount;
            c.hydrogenCount = this->hydrogenCount + a.hydrogenCount;
            c.nitrogenCount = this->nitrogenCount + a.nitrogenCount;
            c.oxygenCount = this->oxygenCount + a.oxygenCount;
            c.sulfurCount = this->sulfurCount + a.sulfurCount;
            c.phosphorusCount = this->phosphorusCount + a.phosphorusCount;

            return c;
        }

        MolecularFormula& operator+=(const MolecularFormula &a) {
            this->carbonCount += a.carbonCount;
            this->hydrogenCount += a.hydrogenCount;
            this->nitrogenCount += a.nitrogenCount;
            this->oxygenCount += a.oxygenCount;
            this->sulfurCount += a.sulfurCount;
            this->phosphorusCount += a.phosphorusCount;

            return *this;
        }

        MolecularFormula& operator-=(const MolecularFormula &a) {
            this->carbonCount -= a.carbonCount;
            this->hydrogenCount -= a.hydrogenCount;
            this->nitrogenCount -= a.nitrogenCount;
            this->oxygenCount -= a.oxygenCount;
            this->sulfurCount -= a.sulfurCount;
            this->phosphorusCount -= a.phosphorusCount;

            return *this;
        }

    };

    Err CHEMLIB_EXPORTS parseMolecularFormulaString(const QString &formulaString, MolecularFormula *mf);

    //Amino Acids
    extern const MolecularFormula CHEMLIB_EXPORTS alanineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS cysteineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS asparticAcidFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS glutamicAcidFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS phenylalanineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS glycineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS histidineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS isoleucineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS lysineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS leucineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS methionineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS asparagineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS prolineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS glutamineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS arginineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS serineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS threonineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS valineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS tryptophanFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS tyrosineFormula;

    //Small Groups CHEMLIB_EXPORTS
    extern const MolecularFormula CHEMLIB_EXPORTS ammoniaFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS carbonylFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS carbamidomethylFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS carboxymethylFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS methyenelFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS methylFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS phosphateFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS waterFormula;

    //Bases CHEMLIB_EXPORTS
    extern const MolecularFormula CHEMLIB_EXPORTS adenineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS cytosineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS guanineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS thymineFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS uracilFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS hypoxanthine;

    //Sugars CHEMLIB_EXPORTS
    extern const MolecularFormula CHEMLIB_EXPORTS deOxyRiboseFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS riboseFormula;

    //Modifications CHEMLIB_EXPORTS
    extern const MolecularFormula CHEMLIB_EXPORTS hydrogenLossFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS oxidationFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS doubleOxidationFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS deamidationFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS ammoniaLossFormula;
    extern const MolecularFormula CHEMLIB_EXPORTS waterLossFormula;

}//NAMESPACE

#endif // MOLECULARFORMULAS_H
