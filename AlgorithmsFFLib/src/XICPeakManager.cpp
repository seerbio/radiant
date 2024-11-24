//
// Created by anichols on 5/9/24.
//

#include "XICPeakManager.h"

#include "EigenKernelUtils.h"
#include "ErrorUtils.h"
#include "MS2Ion.h"
#include "MsFrame.h"
#include "TargetDecoyCandidatePair.h"


XICPeakManager::XICPeakManager()
: m_filterLength(-1)
, m_polynomialOrder(-1)
, m_smoothCount(-1)
, m_stopThresholdValue(-1.0)
, m_isInit(false)
, m_turboXic(nullptr)
, m_deleteTurboXicLocal(false)
, m_ppmTolerance(-1.0)
{}

XICPeakManager::~XICPeakManager() {
    if (m_deleteTurboXicLocal) {
        delete m_turboXic;
    }
}

Err XICPeakManager::init(
    const MsFrame &msFrame,
    const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
    int topNMs2Ions,
    float ppmTolerance
    ){

    ERR_INIT

    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(targetDecoyPointers); ree;
    e = ErrorUtils::isAboveThreshold(ppmTolerance, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    m_ppmTolerance = ppmTolerance;

    m_turboXic = new TurboXIC();
    e = m_turboXic->init(msFrame.frameIndexVsScanPoints()); ree;

    m_isInit = true;
    m_deleteTurboXicLocal = true;

    ERR_RETURN
}

Err XICPeakManager::init(
    const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
    int topNMs2Ions,
    float ppmTolerance,
    TurboXIC* turboXic
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(turboXic->isInit()); ree;
    e = ErrorUtils::isNotEmpty(targetDecoyPointers); ree;
    e = ErrorUtils::isAboveThreshold(ppmTolerance, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    m_ppmTolerance = ppmTolerance;

    m_turboXic = turboXic;

    m_isInit = true;

    ERR_RETURN
}

bool XICPeakManager::isValid() const {
    return m_isInit;
}

Err XICPeakManager:: getXIC(
    float mzVal,
    XICPoints *xicPoints
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isAboveThreshold(mzVal, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    const float massTol = MathUtils::calculatePPM(mzVal, m_ppmTolerance);
    const float mzMin = mzVal - massTol;
    const float mzMax = mzVal + massTol;

    const XICPoints xicPointsLocal = m_turboXic->extractPointsXIC(mzMin, mzMax);
    *xicPoints = xicPointsLocal;

    ERR_RETURN
}

