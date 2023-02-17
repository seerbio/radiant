#ifndef CHEMCONSTANTS_H
#define CHEMCONSTANTS_H


namespace ChemConstants
{

    namespace MonoIsotopicMass {
        extern const double CARBON;
        extern const double HYDROGEN;
        extern const double NITROGEN;
        extern const double OXYGEN;
        extern const double PHOSPHORUS;
        extern const double SULFUR;
    }//END NAMESPACE

    namespace AverageMass {
        extern const double CARBON;
        extern const double HYDROGEN;
        extern const double NITROGEN;
        extern const double OXYGEN;
        extern const double PHOSPHORUS;
        extern const double SULFUR;
    }//END NAMESPACE

    extern const double	PROTON;
    extern const double	NEUTRON;

    extern const double C12_ABUNDANCE;
    extern const double C13_ABUNDANCE;

    extern const double H1_ABUNDANCE;
    extern const double H2_ABUNDANCE;

    extern const double N14_ABUNDANCE;
    extern const double N15_ABUNDANCE;

    extern const double O16_ABUNDANCE;
    extern const double O17_ABUNDANCE;
    extern const double O18_ABUNDANCE;

    extern const double S32_ABUNDANCE;
    extern const double S33_ABUNDANCE;
    extern const double S34_ABUNDANCE;
    extern const double S36_ABUNDANCE;

    struct AVERAGINE {

        const double C = 4.94;
        const double H = 7.76;
        const double N = 1.36;
        const double O = 1.48;
        const double S = 0.04;

        const double averagineMass = C * AverageMass::CARBON
                + H * AverageMass::HYDROGEN
                + N * AverageMass::NITROGEN
                + O * AverageMass::OXYGEN
                + S * AverageMass::SULFUR;
    };

}//NAMESPACE

#endif // CHEMCONSTANTS_H
