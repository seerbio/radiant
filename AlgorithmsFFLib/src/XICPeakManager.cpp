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

namespace {

    Err buildMzHashedVsCount(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
        int topNFragIons,
        QHash<MzHashed, Occurrence> *mzHashedVsCount
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyPointers); eee_absorb;
        for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {

            int counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsTarget()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(
                    ms2Ion.mz,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );
                (*mzHashedVsCount)[mzHashed]++;

                const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / static_cast<float>(ms2Ion.charge);
                const MzHashed mzHashedShadow = MathUtils::hashDecimal(
                    ms2Ion.mz - isotopeDistanceThomsons,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );
                (*mzHashedVsCount)[mzHashedShadow]++;
            }

            counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsDecoy()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;

                const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / static_cast<float>(ms2Ion.charge);
                const MzHashed mzHashedShadow = MathUtils::hashDecimal(
                    ms2Ion.mz - isotopeDistanceThomsons,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );

                (*mzHashedVsCount)[mzHashedShadow]++;

            }
        }

        ERR_RETURN
    }
}//namespace
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

    e = buildMzHashedVsCount(
        targetDecoyPointers,
        topNMs2Ions,
        &m_mzHashedOccurrences
        ); ree;
    e = ErrorUtils::isNotEmpty(m_mzHashedOccurrences); ree;

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

    e = buildMzHashedVsCount(
        targetDecoyPointers,
        topNMs2Ions,
        &m_mzHashedOccurrences
        ); ree;
    e = ErrorUtils::isNotEmpty(m_mzHashedOccurrences); ree;

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
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isAboveThreshold(mzVal, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    // const MzHashed mzHashed = MathUtils::hashDecimal(mzVal, S_GLOBAL_SETTINGS.HASHING_PRECISION);
    //
    // if (m_mzHashedVsXicPoints.contains(mzHashed)) {
    //
    //     *xicPoints = m_mzHashedVsXicPoints.value(mzHashed);
    //
    //     m_mzHashedOccurrences[mzHashed]--;
    //     if (m_mzHashedOccurrences[mzHashed] == 0) {
    //         XICPoints().swap(m_mzHashedVsXicPoints[mzHashed]);
    //         m_mzHashedVsXicPoints.remove(mzHashed);
    //     }
    //     ERR_RETURN
    // }

    const float massTol = MathUtils::calculatePPM(mzVal, m_ppmTolerance);
    const float mzMin = mzVal - massTol;
    const float mzMax = mzVal + massTol;

    const XICPoints xicPointsLocal = m_turboXic->extractPointsXIC(mzMin, mzMax);
    *xicPoints = xicPointsLocal;

    // if (constexpr int minOccuranceCountToCache = 2; m_mzHashedOccurrences.value(mzHashed) < minOccuranceCountToCache) {
    //     ERR_RETURN
    // }
    //
    // m_mzHashedVsXicPoints.insert(mzHashed , xicPointsLocal);
    // m_mzHashedOccurrences[mzHashed]--;

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
