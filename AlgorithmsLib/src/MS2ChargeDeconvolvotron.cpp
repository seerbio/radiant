//
// Created by anichols on 7/7/23.
//

#include "MS2ChargeDeconvolvotron.h"

#include "ErrorUtils.h"
#include "NeuralNetModel.h"


class Q_DECL_HIDDEN MS2ChargeDeconvolvotron::Private
{
public:

    Private();
    ~Private();

    Err init(
            const QString &modelFilePath,
            double ppmTolerance
    );

    Err predictMS2Charge(
            const ScanPoints &scanPoints,
            int *charge
            );

private:

    NeuralNetModel* m_model;
    int m_chargeMin;
    int m_chargeMax;
    double m_ppmTolerance;
    bool m_isInit;

};

MS2ChargeDeconvolvotron::Private::Private()
: m_chargeMin(1)
, m_chargeMax(3)
, m_ppmTolerance(10.0)
, m_model(nullptr)
, m_isInit(false)
{}

MS2ChargeDeconvolvotron::Private::~Private() {
    if(m_isInit){
        delete m_model;
    }
}

Err MS2ChargeDeconvolvotron::Private::init(
        const QString &modelFilePath,
        double ppmTolerance
        ) {

    ERR_INIT

    const double ppmToleranceMin = 0.0;

    e = ErrorUtils::fileExists(modelFilePath); ree
    e = ErrorUtils::isAboveThreshold(
            ppmTolerance,
            ppmToleranceMin,
            ErrorUtilsParam::ExcludeThreshold
            ); ree

    m_ppmTolerance = ppmTolerance;

    m_model = new NeuralNetModel();
    e = m_model->init(modelFilePath); ree;

    m_isInit = true;

    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::Private::predictMS2Charge(
        const ScanPoints &scanPoints,
        int *charge
        ) {
    return eFunctionNotImplemented;
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


MS2ChargeDeconvolvotron::MS2ChargeDeconvolvotron()
: d_ptr(new Private()) {}

MS2ChargeDeconvolvotron::~MS2ChargeDeconvolvotron() {}

Err MS2ChargeDeconvolvotron::init(
        const QString &modelFilePath,
        double ppmTolerance
        ) {
    ERR_INIT
    e = d_ptr->init(modelFilePath, ppmTolerance); ree
    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::predictMS2Charge(
        const ScanPoints &scanPoints,
        int *charge
        ) {
    ERR_INIT
    e = d_ptr->predictMS2Charge(scanPoints, charge); ree
    ERR_RETURN
}
