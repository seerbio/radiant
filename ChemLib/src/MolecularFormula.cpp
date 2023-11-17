#include "MolecularFormula.h"

#include "ErrorUtils.h"

#include <QDebug>
#include <QRegExp>
#include <QString>
#include <QStringList>

namespace MolecularFormulas
{

const MolecularFormula alanineFormula(3, 5, 1, 1, 0, 0); //C3H5NO1
const MolecularFormula cysteineFormula(3, 5, 1, 1, 1, 0); //C3H5NOS
const MolecularFormula asparticAcidFormula(4, 5, 1, 3, 0, 0); //C4H5NO3
const MolecularFormula glutamicAcidFormula(5, 7, 1, 3, 0, 0); //C5H7NO3
const MolecularFormula phenylalanineFormula(9, 9, 1, 1, 0, 0); //C9H9NO
const MolecularFormula glycineFormula(2, 3, 1, 1, 0, 0); //C2H3NO
const MolecularFormula histidineFormula(6, 7, 3, 1, 0, 0); //C6H7N3O
const MolecularFormula isoleucineFormula(6, 11, 1, 1, 0, 0); //C6H11NO
const MolecularFormula lysineFormula(6, 12, 2, 1, 0, 0); //C6H12N2O
const MolecularFormula leucineFormula(6, 11, 1, 1, 0, 0); //C6H11NO
const MolecularFormula methionineFormula(5, 9, 1, 1, 1, 0); //C5H9NOS
const MolecularFormula asparagineFormula(4, 6, 2, 2, 0, 0); //C4H6N2O2
const MolecularFormula prolineFormula(5, 7, 1, 1, 0, 0); //C5H7NO
const MolecularFormula glutamineFormula(5, 8, 2, 2, 0, 0); //C5H8N2O2
const MolecularFormula arginineFormula(6, 12, 4, 1, 0, 0); //C6H12N4O
const MolecularFormula serineFormula(3, 5, 1, 2, 0, 0); //C3H5NO2
const MolecularFormula threonineFormula(4, 7, 1, 2, 0, 0); //C4H7NO2
const MolecularFormula valineFormula(5, 9, 1, 1, 0, 0); //C5H9NO
const MolecularFormula tryptophanFormula(11, 10, 2, 1, 0, 0); //C11H10N2O
const MolecularFormula tyrosineFormula(9, 9, 1, 2, 0, 0); //C9H9NO2

const MolecularFormula ammoniaFormula(0, 3, 1, 0, 0, 0);//H3N
const MolecularFormula carbonylFormula(1, 0, 0, 1, 0, 0);//CO
const MolecularFormula carbamidomethylFormula(2, 3, 1, 1, 0, 0); //C2H3N1O1
const MolecularFormula carboxymethylFormula(2, 2, 0, 2, 0, 0); //C2H2NO2;
const MolecularFormula methyenelFormula(2, 4, 0, 0, 0, 0);//C2H4;
const MolecularFormula methylFormula(1, 3, 0, 0, 0, 0);//CH3
const MolecularFormula phosphateFormula(0, 0, 0, 4, 0, 1);//O4P;
const MolecularFormula waterFormula(0, 2, 0, 1, 0, 0); //H2O

const MolecularFormula adenineFormula(5, 5, 5, 0, 0, 0);//C5H5N5
const MolecularFormula cytosineFormula(4, 5, 3, 1, 0, 0); //C4H5N3O
const MolecularFormula guanineFormula(5, 5, 5, 1, 0, 0);//C5H5N5O
const MolecularFormula thymineFormula(5, 6, 2, 2, 0, 0);//C5H6N2O2
const MolecularFormula uracilFormula(4, 4, 2, 2, 0, 0);//C4H4N2O2
const MolecularFormula hypoxanthine(5, 4, 4, 1, 0, 0);//C5H4N4O

const MolecularFormula deOxyRiboseFormula(5, 10, 0, 4, 0, 0);//C5H10O4
const MolecularFormula riboseFormula(5, 10, 0, 5, 0, 0);//C5H10O5

/***********************************
 * Modifications
 **********************************/
//CHNOS
const MolecularFormula hydrogenLossFormula(0, -1, 0, 0, 0, 0);
const MolecularFormula oxidationFormula(0, 0, 0, 1, 0, 0);
const MolecularFormula doubleOxidationFormula(0, 0, 0, 2, 0, 0);
const MolecularFormula deamidationFormula(0, -1, -1, 1, 0, 0);
const MolecularFormula ammoniaLossFormula(0, -3, -1, 0, 0, 0);
const MolecularFormula waterLossFormula(0, -2, 0, -1, 0, 0);

Err CHEMLIB_EXPORTS parseMolecularFormulaString(const QString &_formulaString, MolecularFormula *mf) {

    ERR_INIT

    const QString allowedAtoms = QStringLiteral("CHNOSP");
    QString formulaString = _formulaString.toUpper();
    formulaString = formulaString.remove("(");
    formulaString = formulaString.remove(")");
    formulaString = formulaString.remove(" ");

    QHash<QString, int> formulaHash;
    for (const QString &c : allowedAtoms.split("", Qt::SkipEmptyParts)) {
        formulaHash.insert(c, 0);
    }

    QStringList molStringSplit = formulaString.split("", Qt::SkipEmptyParts);

    const QRegExp rx(QLatin1String("[-0-9]+"));
    auto&& atomsInString
            = formulaString.split(rx, Qt::SkipEmptyParts).join("").split("", Qt::SkipEmptyParts);

    const auto logic = [allowedAtoms](const QString &c){return !allowedAtoms.contains(c); };
    const bool areIllegalAtomsInFormula = std::any_of(atomsInString.begin(), atomsInString.end(), logic);
    e = ErrorUtils::isFalse(areIllegalAtomsInFormula); ree;

    QString currentAtom;
    QString numberBuilder;
    while (!molStringSplit.isEmpty()) {

        const QString &currentFront =  molStringSplit.front();
        if (atomsInString.contains(currentFront)) {

            if (!currentAtom.isEmpty()) {
                bool ok = true;
                const int atomCount = numberBuilder.isEmpty() ? 1 : numberBuilder.toInt(&ok);
                e = ErrorUtils::isTrue(ok); ree;

                formulaHash[currentAtom] = atomCount;
                numberBuilder.clear();
            }

            currentAtom = currentFront;
            atomsInString.pop_front();
            molStringSplit.pop_front();
            continue;
        }

        numberBuilder += currentFront;
        molStringSplit.pop_front();
    }

    bool ok = true;
    const int atomCount = numberBuilder.isEmpty() ? 1 : numberBuilder.toInt(&ok);
    e = ErrorUtils::isTrue(ok); ree;
    formulaHash[currentAtom] = atomCount;

    *mf = MolecularFormula(formulaHash["C"],
                           formulaHash["H"],
                           formulaHash["N"],
                           formulaHash["O"],
                           formulaHash["S"],
                           formulaHash["P"]);

    ERR_RETURN
}


}//NAMESPACE
