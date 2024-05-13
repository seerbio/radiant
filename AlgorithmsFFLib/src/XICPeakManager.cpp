//
// Created by anichols on 5/9/24.
//

#include "XICPeakManager.h"

#include "EigenKernelUtils.h"
#include "ErrorUtils.h"
#include "MsFrame.h"


XICPeakManager::XICPeakManager()
: m_filterLength(-1)
, m_polynomialOrder(-1)
, m_smoothCount(-1)
, m_stopThresholdValue(-1.0)
{}

Err XICPeakManager::init(
    const MsFrame &msFrame,
    const QVector<float> &mzValsToExtract,
    float ppmTolerance
    ){

    ERR_INIT

    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(mzValsToExtract); ree;
    e = ErrorUtils::isAboveThreshold(ppmTolerance, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    TurboXIC turboXic;
    e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

    e = extractXICs(
        mzValsToExtract,
        ppmTolerance,
        &turboXic
        ); ree;

    ERR_RETURN
}

bool XICPeakManager::isValid() const {
    const QList<MzHashed> &keys = m_mzHashedVsXicPoints.keys();
    return !m_mzHashedVsXicPoints.isEmpty() && !m_mzHashedVsXicPoints.value(keys.first()).empty();
}

Err XICPeakManager:: getXIC(
    float mzVal,
    XICPoints *xicPoints
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_mzHashedVsXicPoints); ree;
    e = ErrorUtils::isAboveThreshold(mzVal, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    const MzHashed mzHashed = MathUtils::hashDecimal(mzVal, S_GLOBAL_SETTINGS.HASHING_PRECISION);

    *xicPoints = m_mzHashedVsXicPoints.value(mzHashed);
    ERR_RETURN
}

 Err XICPeakManager::cacheXICPeakManager(const QString &outputFilePath) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_mzHashedVsXicPoints); ree;

    QFile file(outputFilePath);

    if(file.exists()) {
        const bool removed = file.remove();
        e = ErrorUtils::isTrue(removed);
    }

    if (!file.open(QIODevice::WriteOnly)) {
        qDebug() << "Failed to open file!";
        rrr(eFileError);
    }

    QDataStream out(&file);
    out.setVersion(QDataStream::Qt_5_15);

    const QList<XICPoints> &xicPointses = m_mzHashedVsXicPoints.values();

    QVector<QVector<XICPoint>> vec;
    std::transform (
        xicPointses.begin(),
        xicPointses.end(),
        std::back_inserter(vec),
        [](const XICPoints &p){return QVector<XICPoint>(p.begin(), p.end());}
        );

    out << m_filterLength;
    out << m_polynomialOrder;
    out << m_smoothCount;
    out << m_stopThresholdValue;
    out << m_mzHashedVsXicPoints.keys().toVector();
    out << vec;

    file.close();

    ERR_RETURN
}

Err XICPeakManager::loadXICPeakManagerCache(const QString& outputFilePath) {

    ERR_INIT
    e = ErrorUtils::fileExists(outputFilePath); ree;

    QFile file(outputFilePath);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open file!";
        rrr(eFileError);
    }

    QDataStream in(&file);
    in.setVersion(QDataStream::Qt_5_15);

    QVector<MzHashed> mzHashedVec;
    QVector<QVector<XICPoint>> vec;

    in >> m_filterLength;
    in >> m_polynomialOrder;
    in >> m_smoothCount;
    in >> m_stopThresholdValue;
    in >> mzHashedVec;
    in >> vec;

    e = ErrorUtils::isEqual(mzHashedVec.size(), vec.size()); ree;
    for (int i = 0; i < mzHashedVec.size(); i++) {
        const QVector<XICPoint> &v = vec.at(i);
        m_mzHashedVsXicPoints.insert(mzHashedVec.at(i), std::vector<XICPoint>(v.begin(), v.end()));
    }

    file.close();

    ERR_RETURN
}

Err XICPeakManager::extractXICs(
    const QVector<float> &mzValsToExtract,
    float ppmTolerance,
    TurboXIC *turboXic
) {

    ERR_INIT

    for (float mzVal : mzValsToExtract) {

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTolerance);
        const float mzMin = mzVal - massTol;
        const float mzMax = mzVal + massTol;

        const XICPoints xicPoints = turboXic->extractPointsXIC(mzMin, mzMax);

        if (xicPoints.empty()) {
            continue;
        }

        const MzHashed mzHashed = MathUtils::hashDecimal(mzVal, S_GLOBAL_SETTINGS.HASHING_PRECISION);
        m_mzHashedVsXicPoints.insert(mzHashed , xicPoints);

    }

    ERR_RETURN
}
