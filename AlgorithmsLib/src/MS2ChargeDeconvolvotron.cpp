//
// Created by anichols on 7/7/23.
//

#include "MS2ChargeDeconvolvotron.h"

#include "BiophysicalCalcs.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "IsotopicDistributionBuilder.h"
#include "MS2PointsExtractomatic.h"
#include "MsUtils.h"
#include "NeuralNetModel.h"
#include "ParallelUtils.h"


class Q_DECL_HIDDEN MS2ChargeDeconvolvotron::Private
{
public:

    using MassInt = int;

    Private();
    ~Private();

    Err init(
            const QString &chargeModelFilePath,
            const QString &monoModelFilePath,
            double ppmTolerance
    );

    Err buildMzIntVsAveragineRatios();

    Err deisotopeScanPoints(
            const ScanPoints &scanPoints,
            ScanPoints *scanPointsDeisotoped
            );

    QVector<double> buildMzValsForChargeMono(double mzMaxIntensity) const;

    Err extractScanPointsFromScan(
            const QVector<double> &mzValsToExtract,
            MS2PointsExtractomatic *ms2PointsExtractomatic,
            QVector<QPair<Id, ScanPoint>> *foundScanPoints
    ) const;

    Err predictChargeAnMonoOffset(
            const QVector<QPair<Id, ScanPoint>> &foundScanPoints,
            int *charge,
            int *monoOffset,
            QVector<float> *predVec,
            QVector<float> *predVecMono
    );

    Err testChargeMonoCaller(
            const ScanPoints &scanPoints,
            double mzVal,
            QVector<float> *predVecCharge,
            QVector<float> *predVecMono
            );

    Err buildSubtractionPoints(
            const QVector<QPair<Id, ScanPoint>> &chargePoints,
            int charge,
            QVector<QPair<Id, ScanPoint>> *chargePointsSubtracted
            );

private:

    NeuralNetModel *m_modelCharge;
    NeuralNetModel *m_modelMonoOffset;
    int m_chargeMin;
    int m_chargeMax;
    double m_ppmTolerance;
    bool m_isInit;
    QMap<MassInt , Eigen::VectorX<double>> m_mzIntVsAveragineRatios;

};

MS2ChargeDeconvolvotron::Private::Private()
: m_chargeMin(1)
, m_chargeMax(3)
, m_ppmTolerance(15.0)
, m_modelCharge(nullptr)
, m_modelMonoOffset(nullptr)
, m_isInit(false)
{}

MS2ChargeDeconvolvotron::Private::~Private() {
    if(m_isInit){
        delete m_modelCharge;
        delete m_modelMonoOffset;
    }
}

Err MS2ChargeDeconvolvotron::Private::init(
        const QString &chargeModelFilePath,
        const QString &monoModelFilePath,
        double ppmTolerance
        ) {

    ERR_INIT

    const double ppmToleranceMin = 0.0;

    e = ErrorUtils::fileExists(chargeModelFilePath); ree
    e = ErrorUtils::fileExists(monoModelFilePath); ree
    e = ErrorUtils::isAboveThreshold(
            ppmTolerance,
            ppmToleranceMin,
            ErrorUtilsParam::ExcludeThreshold
            ); ree

    m_ppmTolerance = ppmTolerance;

    m_modelCharge = new NeuralNetModel();
    e = m_modelCharge->init(chargeModelFilePath); ree

    m_modelMonoOffset = new NeuralNetModel();
    e = m_modelMonoOffset->init(monoModelFilePath); ree

    e = buildMzIntVsAveragineRatios(); ree

    m_isInit = true;

    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::Private::buildMzIntVsAveragineRatios() {

    ERR_INIT

    const int massMin = 50;
    const int massMax = 5000;

    for (MassInt mass = massMin; mass < massMax; mass++) {

        const QVector<double> isoDistAtMass
                = IsotopicDistributionBuilder::getIsotopicDistribution(static_cast<double>(mass));

        const Eigen::VectorX<double> isoDistAtMassVec
                = EigenUtils::convertQVectorToEigenVector(isoDistAtMass);

        m_mzIntVsAveragineRatios.insert(mass, isoDistAtMassVec);
    }

    e = ErrorUtils::isNotEmpty(m_mzIntVsAveragineRatios); ree

    ERR_RETURN
}

namespace {

    const auto maxLogicIntensity = [](const QPair<Id, ScanPoint> &l, const QPair<Id, ScanPoint> &r){
        return l.second.y() < r.second.y();
    };

    void sortIdScanPointPairsHiLoIntensity(QVector<QPair<Id, ScanPoint>> *idScanPointPairsHiLoIntensity) {

        std::sort(idScanPointPairsHiLoIntensity->rbegin(), idScanPointPairsHiLoIntensity->rend(), maxLogicIntensity);
    }

    QVector<QPair<Id, ScanPoint>> idVsScanPointsToPairsHiLoIntensity(const QMap<Id, ScanPoint> &idVsScanPoints) {

        QVector<QPair<Id, ScanPoint>> idScanPointPairsHiLoIntensity;
        for (auto it = idVsScanPoints.begin(); it != idVsScanPoints.end(); it++) {

            const Id id = it.key();
            const ScanPoint &sp = it.value();

            idScanPointPairsHiLoIntensity.append({id, sp});
        }

        sortIdScanPointPairsHiLoIntensity(&idScanPointPairsHiLoIntensity);

        return idScanPointPairsHiLoIntensity;
    }

    Err extractScanPointsByCharge(
            const QVector<QPair<Id, ScanPoint>> &extractedPoints,
            int charge,
            int monoOffset,
            QVector<QPair<Id, ScanPoint>> *outScanPoints
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(extractedPoints); ree

        int idxStart = -1;
        int idxEnd = -1;

        switch (charge) {

            case 1:
                idxStart = 0 + charge - monoOffset;
                idxEnd = 4;
                break;
            case 2:
                idxStart = 4 + charge - monoOffset;
                idxEnd = 10;
                break;
            case 3:
                idxStart = 10 + charge - monoOffset;
                idxEnd = 18;
                break;
            default:
                idxStart = 0 + charge - monoOffset;
                idxEnd = 4;
        }

        *outScanPoints = extractedPoints.mid(idxStart, idxEnd - idxStart);

        ERR_RETURN
    }

    bool isMaxIntensityZero(const QVector<QPair<Id, ScanPoint>> &chargePoints) {

        const double chargePointsMaxIntensity = std::max_element(
                chargePoints.begin(),
                chargePoints.end(),
                maxLogicIntensity
        )->second.y();

        return MathUtils::tZero(chargePointsMaxIntensity) || chargePointsMaxIntensity < 0;

    }

    const auto extractPointsTransformLogic = [](const QPair<Id, ScanPoint> &p){
        return p.second;
    };

    Err updateScanPoints(
            const QVector<QPair<Id, ScanPoint>> &chargePointsSubtracted,
            MS2PointsExtractomatic *ms2PointsExtractomatic
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(chargePointsSubtracted); ree

        QMap<Id, ScanPoint> scanPointsMapUpdate;
        QMap<Id, ScanPoint> scanPointsMapDelete;
        for (const QPair<Id, ScanPoint> &sp : chargePointsSubtracted) {

            const Id id = sp.first;
            const ScanPoint &p = sp.second;

            if (id < 0) {
                continue;
            }

            if (p.y() <= 0) {
                scanPointsMapDelete.insert(id, p);
                continue;
            }

            scanPointsMapUpdate.insert(id, p);
        }

        if (!scanPointsMapUpdate.isEmpty()) {
            e = ms2PointsExtractomatic->updateScanPoints(scanPointsMapUpdate); ree
        }

        if (!scanPointsMapDelete.isEmpty()) {
            e = ms2PointsExtractomatic->removeScanPoints(scanPointsMapDelete); ree
        }

        ERR_RETURN
    }

}//namespace
Err MS2ChargeDeconvolvotron::Private::deisotopeScanPoints(
        const ScanPoints &scanPoints,
        ScanPoints *scanPointsDeisotoped
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init({scanPoints}); ree

    const QMap<Id, ScanPoint> idVsScanPoints = ms2PointsExtractomatic.idVsScanPoints();

    const QVector<QPair<Id, ScanPoint>> idScanPointPairsHiLoIntensity = idVsScanPointsToPairsHiLoIntensity(idVsScanPoints);

    for (const QPair<Id, ScanPoint> &idsp : idScanPointPairsHiLoIntensity) {

        const double mzExtract = idsp.second.x();

        const QVector<double> mzValsToExtract = buildMzValsForChargeMono(mzExtract);

        QVector<QPair<Id, ScanPoint>> extractedMzExtract;
        e = extractScanPointsFromScan(
                {mzExtract},
                &ms2PointsExtractomatic,
                &extractedMzExtract
        ); ree

        const QPair<Id, ScanPoint> &mzExtractScanPoint = extractedMzExtract.front();
        if (mzExtractScanPoint.second.y() < 0 || MathUtils::tZero(mzExtractScanPoint.second.y())) {
            continue;
        }

        QVector<QPair<Id, ScanPoint>> extractedScanPoints;
        e = extractScanPointsFromScan(
                mzValsToExtract,
                &ms2PointsExtractomatic,
                &extractedScanPoints
                ); ree

        e = ErrorUtils::isEqual(mzValsToExtract.size(), extractedScanPoints.size()); ree

        int charge;
        int monoOffset;
        QVector<float> unusedPredVecCharge;
        QVector<float> unusedPredVecMono;
        e = predictChargeAnMonoOffset(
                extractedScanPoints,
                &charge,
                &monoOffset,
                &unusedPredVecCharge,
                &unusedPredVecMono
                ); ree

        qDebug() << mzExtract << charge << monoOffset;

        if (charge < 1 || charge == 3) {
            continue;
        }

        if (charge == 1) { //TODO, fix the model so you don't need this.
            monoOffset = 0;
        }

        QVector<QPair<Id, ScanPoint>> chargePoints;
        e = extractScanPointsByCharge(
                extractedScanPoints,
                charge,
                monoOffset,
                &chargePoints
                ); ree

        if (isMaxIntensityZero(chargePoints)) {
            continue;
        }

        const QPair<Id, ScanPoint> &chargePointsMonoIsotope = chargePoints.at(0);
        const double mzMin = 1.0;
        if (chargePointsMonoIsotope.second.x() < mzMin || MathUtils::tZero(chargePointsMonoIsotope.second.y())) {
            continue;
        }

        scanPointsDeisotoped->push_back(chargePointsMonoIsotope.second);

        QVector<QPair<Id, ScanPoint>> chargePointsSubtracted;
        e = buildSubtractionPoints(chargePoints, charge, &chargePointsSubtracted); ree
        e = updateScanPoints(chargePointsSubtracted, &ms2PointsExtractomatic); ree
    }

    std::sort(scanPointsDeisotoped->begin(), scanPointsDeisotoped->end(), MsUtilsNamespace::sortAscMz);

    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::Private::buildSubtractionPoints(
        const QVector<QPair<Id, ScanPoint>> &chargePoints,
        int charge,
        QVector<QPair<Id, ScanPoint>> *chargePointsSubtracted
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(chargePoints); ree

    chargePointsSubtracted->clear();

    const double subtractionBufferFactor = 1.2;

    ScanPoints scanPoints;
    std::transform(chargePoints.begin(), chargePoints.end(), std::back_inserter(scanPoints), extractPointsTransformLogic);

    QPair<QVector<double>, QVector<double>> mzValsintensityValsSplit = ParallelUtils::unZip(scanPoints);
    const QVector<double> &mzVals = mzValsintensityValsSplit.first;
    const QVector<double> &intensityVals = mzValsintensityValsSplit.second;
    const double intensityValsMax = *std::max_element(intensityVals.begin(), intensityVals.end());


    const auto roughMass = static_cast<MassInt>(mzVals.at(0) * charge);
    Eigen::VectorX<double> mzAveragine = m_mzIntVsAveragineRatios.value(roughMass);
    mzAveragine.conservativeResize(chargePoints.size());
    mzAveragine *= (subtractionBufferFactor * intensityValsMax);

    Eigen::VectorX<double> intensityValsVec = EigenUtils::convertQVectorToEigenVector(intensityVals);
    intensityValsVec -= mzAveragine;

    e = ErrorUtils::isEqual(chargePoints.size(), intensityVals.size()); ree

    for (int i = 0; i < chargePoints.size(); i++) {

        const double intensity = intensityValsVec.coeff(i);
        QPair<Id, ScanPoint> subtractedPoint = chargePoints.at(i);

        subtractedPoint.second.ry() = intensity;

        chargePointsSubtracted->push_back(subtractedPoint);
    }

    ERR_RETURN
}

QVector<double> MS2ChargeDeconvolvotron::Private::buildMzValsForChargeMono(double mzMaxIntensity) const {

    QVector<double> extractMzVals;

    for (int charge = m_chargeMin; charge <= m_chargeMax; charge++) {

        const QVector<double> extractMzValsCharge = BiophysicalCalcs::calculateIsotopesFromMz(
                mzMaxIntensity,
                charge,
                charge,
                charge + 1
                );

        extractMzVals.append(extractMzValsCharge);
    }

    return extractMzVals;
}

namespace {

    QPair<Id, ScanPoint> returnClosestScanPointToMz(
            const QMap<Id, ScanPoint> &idVsScanPoint,
            double mz
            ) {

        QPair<Id, ScanPoint> bestScanPoint = {-1, {std::numeric_limits<double>::max() , -1}};

        for (auto it = idVsScanPoint.begin(); it != idVsScanPoint.end(); it++) {

            const Id id = it.key();
            const ScanPoint &sp = it.value();

            const double mzDiffBestScanPoint = std::abs(mz - bestScanPoint.second.x());
            const double mzDiffCurrentPoint = std::abs(mz - sp.x());

            if (mzDiffCurrentPoint < mzDiffBestScanPoint) {
                bestScanPoint = {id, sp};
            }

        }

        return bestScanPoint;
    }

}//namespace
Err MS2ChargeDeconvolvotron::Private::extractScanPointsFromScan(
        const QVector<double> &mzValsToExtract,
        MS2PointsExtractomatic *ms2PointsExtractomatic,
        QVector<QPair<Id, ScanPoint>> *foundScanPoints
) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzValsToExtract); ree

    for (double mz : mzValsToExtract) {

        const double mzTolerance = MathUtils::calculatePPM(mz, m_ppmTolerance);
        const double mzMin = mz - mzTolerance;
        const double mzMax = mz + mzTolerance;

        QMap<Id, ScanPoint> foundScanPointsMz;
        e = ms2PointsExtractomatic->findScanPoints(
                mzMin,
                mzMax,
                &foundScanPointsMz
                ); ree

        if (foundScanPointsMz.isEmpty()) {
            foundScanPoints->push_back({-1, {mz, 0.0}});
            continue;
        }

        const QPair<Id, ScanPoint> closestMzFoundScanPointsMz = returnClosestScanPointToMz(
                foundScanPointsMz,
                mz
                );

        foundScanPoints->push_back({closestMzFoundScanPointsMz.first, closestMzFoundScanPointsMz.second});

    }

    ERR_RETURN
}

//TODO consider making this into its own class
Err MS2ChargeDeconvolvotron::Private::predictChargeAnMonoOffset(
        const QVector<QPair<Id, ScanPoint>> &foundScanPoints,
        int *charge,
        int *monoOffset,
        QVector<float> *predVecCharge,
        QVector<float> *predVecMono
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree

    *charge = -1;
    *monoOffset = -1;

    ScanPoints scanPoints;
    std::transform(foundScanPoints.begin(), foundScanPoints.end(), std::back_inserter(scanPoints), extractPointsTransformLogic);

    QPair<QVector<double>, QVector<double>> mzValsintensityValsSplit = ParallelUtils::unZip(scanPoints);
    const QVector<double> &mzVals = mzValsintensityValsSplit.first;
    const QVector<double> &intensityVals = mzValsintensityValsSplit.second;

    const double intensityValsMax = *std::max_element(intensityVals.begin(), intensityVals.end());

    if (MathUtils::tZero(intensityValsMax)) {
        ERR_RETURN
    }

    QVector<double> mzValsNorm;
    const double mzNormalizer = 1000.0;
    const auto mzNormalizerLogic = [mzNormalizer](double mz){return mz / mzNormalizer;};
    std::transform(mzVals.begin(), mzVals.end(), std::back_inserter(mzValsNorm), mzNormalizerLogic);

    QVector<double> intensityValsNorm;
    const auto intensityNormalizerLogic = [intensityValsMax](double intz){return intz / intensityValsMax;};
    std::transform(intensityVals.begin(), intensityVals.end(), std::back_inserter(intensityValsNorm), intensityNormalizerLogic);

    QVector<double> mzIntensityCombined;
    mzIntensityCombined.append(mzValsNorm);
    mzIntensityCombined.append(intensityValsNorm);

    const Eigen::VectorX<double> chargeVec = EigenUtils::convertQVectorToEigenVector(intensityValsNorm);
    const Eigen::VectorX<double> monoVec = EigenUtils::convertQVectorToEigenVector(mzIntensityCombined);

//    std::cout << vec << std::endl;

    *predVecCharge = m_modelCharge->predict(chargeVec);
//    qDebug() << intensityValsNorm;
//    qDebug() << *predVecCharge;

    double bestChargeScore = 0;
    for (int i = 0; i < 3; i++) {

        const double scoreAtIndex = predVecCharge->at(i);
        if (scoreAtIndex > bestChargeScore) {
            *charge = i + 1;
            bestChargeScore = scoreAtIndex;
        }
    }

    *predVecMono = m_modelMonoOffset->predict(monoVec);
//    qDebug() << mzIntensityCombined;
//    qDebug() << *predVecMono;

    double bestMonoOffsetScore = 0;

    for (int i = 0; i < predVecMono->size(); i++) {
        const double scoreAtIndex = predVecMono->at(i);
        if (scoreAtIndex > bestMonoOffsetScore) {
            *monoOffset = i;
            bestMonoOffsetScore = scoreAtIndex;
        }
    }

    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::Private::testChargeMonoCaller(
        const ScanPoints &scanPoints,
        double mzVal,
        QVector<float> *predVecCharge,
        QVector<float> *predVecMono
        ) {

    ERR_INIT

    MS2PointsExtractomatic ms2PointsExtractomatic;
    e = ms2PointsExtractomatic.init({scanPoints}); ree

    const QVector<double> mzValsToExtract = buildMzValsForChargeMono(mzVal);

    QVector<QPair<Id, ScanPoint>> foundScanPoints;
    e = extractScanPointsFromScan(
            mzValsToExtract,
            &ms2PointsExtractomatic,
            &foundScanPoints
    ); ree

    e = ErrorUtils::isEqual(mzValsToExtract.size(), foundScanPoints.size()); ree

    int charge;
    int monoOffset;
    e = predictChargeAnMonoOffset(
            foundScanPoints,
            &charge,
            &monoOffset,
            predVecCharge,
            predVecMono
    ); ree

    ERR_RETURN
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


MS2ChargeDeconvolvotron::MS2ChargeDeconvolvotron()
: d_ptr(new Private()) {}

MS2ChargeDeconvolvotron::~MS2ChargeDeconvolvotron() {}

Err MS2ChargeDeconvolvotron::init(
        const QString &chargeModelFilePath,
        const QString &monoModelFilePath,
        double ppmTolerance
        ) {
    ERR_INIT
    e = d_ptr->init(chargeModelFilePath, monoModelFilePath, ppmTolerance); ree
    ERR_RETURN
}

Err MS2ChargeDeconvolvotron::deisotopeScanPoints(
        const ScanPoints &scanPoints,
        ScanPoints *scanPointsDeisotoped
        ) {
    ERR_INIT
    e = d_ptr->deisotopeScanPoints(
            scanPoints,
            scanPointsDeisotoped
            ); ree
    ERR_RETURN
}

QVector<double> MS2ChargeDeconvolvotron::buildMzValsForChargeMono(double mzMaxIntensity) {
    return d_ptr->buildMzValsForChargeMono(mzMaxIntensity);
}

Err MS2ChargeDeconvolvotron::testChargeMonoCaller(
        const ScanPoints &scanPoints,
        double mzVal,
        QVector<float> *predVecCharge,
        QVector<float> *predVecMono
        ) {
    ERR_INIT

    e = d_ptr->testChargeMonoCaller(
            scanPoints,
            mzVal,
            predVecCharge,
            predVecMono
            ); ree

    ERR_RETURN
}
