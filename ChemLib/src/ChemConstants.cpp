#include "ChemConstants.h"
#include "GlobalSettings.h"

namespace ChemConstants
{

    namespace MonoIsotopicMass {
        const double CARBON = 12.0000;
        const double HYDROGEN = 1.007825;
        const double NITROGEN = 14.003074;
        const double OXYGEN = 15.994915;
        const double PHOSPHORUS = 30.973761;
        const double SULFUR = 31.972071;
    }//END NAMESPACE

    namespace AverageMass {
        const double CARBON = 12.0107;
        const double HYDROGEN = 1.00794;
        const double NITROGEN = 14.0067;
        const double OXYGEN = 15.9994;
        const double PHOSPHORUS = 30.9738;
        const double SULFUR = 32.065;
    }//END NAMESPACE

    const double PROTON = 1.00727647;
    const double NEUTRON = S_GLOBAL_SETTINGS.ISO_DIFF;

    const double percentToDecimal = 100.0;

    const double C13_ABUNDANCE = 1.1078 / percentToDecimal;
    const double C12_ABUNDANCE = (1 - C13_ABUNDANCE);

    const double H2_ABUNDANCE = 0.015574 / percentToDecimal;
    const double H1_ABUNDANCE = (1 - H2_ABUNDANCE);

    const double N15_ABUNDANCE = 0.3663 / percentToDecimal;
    const double N14_ABUNDANCE = (1 - N15_ABUNDANCE);

    const double O17_ABUNDANCE = 0.0372 / percentToDecimal;
    const double O18_ABUNDANCE = 0.20004 / percentToDecimal;
    const double O16_ABUNDANCE = (1 - O17_ABUNDANCE - O18_ABUNDANCE);

    const double S33_ABUNDANCE = 0.75 / percentToDecimal;
    const double S34_ABUNDANCE = 4.215 / percentToDecimal;
    const double S36_ABUNDANCE = 0.017 / percentToDecimal;
    const double S32_ABUNDANCE = (1 - S33_ABUNDANCE - S34_ABUNDANCE - S36_ABUNDANCE);

}//NAMESPACE
