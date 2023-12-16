//
// Created by anichols on 10/20/23.
//

#include "CandidateScores.h"

#include "BiophysicalCalcs.h"

Err CandidateScores::buildScoreVector(
        const CandidateScores &candidateScores,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        const QPair<double, double> &scanTimeMinMax,
        QVector<double> *scoreVec
) {

    ERR_INIT

    QVector<double> scores = {
            static_cast<double>(std::max(0, candidateScores.allignedMaxIndexesCount)), //0
            std::max(candidateScores.cosineSimSum100, 0.0), //1
            std::max(candidateScores.cosineSim100MS1, 0.0), //2
            std::pow(std::max(0.0, candidateScores.cosineSimSpectrum), 3), //3
            std::pow(std::max(0.0, candidateScores.klDivSpectrum), 1 / 3.0) //4
    };

//    const QVector<double> columnApexIndexRatiosToAnchorVec
//            = extractScoresFromVecFeatures(candidateScores.columnApexIndexRatiosToAnchor, theoMzIonsSize);
//    scores.append(columnApexIndexRatiosToAnchorVec);

    if (useExtendedScores || useNeuralNetworkScores) {

        const int mzPeakLengthsSum = std::accumulate(
                candidateScores.mzPeakLengthsVec.begin(),
                candidateScores.mzPeakLengthsVec.end(),
                0
        );
        QVector<double> mzPeakLengthsNormalized;
        if (mzPeakLengthsSum != 0) {
            std::transform(
                    candidateScores.mzPeakLengthsVec.begin(),
                    candidateScores.mzPeakLengthsVec.end(),
                    std::back_inserter(mzPeakLengthsNormalized),
                    [mzPeakLengthsSum](int i){return i / static_cast<double>(mzPeakLengthsSum);}
            );
        }
        const QVector<double> mzPeakLengthsNormalizedVec
                = extractScoresFromVecFeatures(mzPeakLengthsNormalized, theoMzIonsSize);
        scores.append(mzPeakLengthsNormalizedVec);

        const double scanTimeDelta = candidateScores.scanTime - candidateScores.scanTimePredicted;
        scores.push_back(std::abs(scanTimeDelta)); //5

        const double scanTimeRange = std::max(
                std::numeric_limits<double>::min(),
                scanTimeMinMax.second - scanTimeMinMax.first
        );

        const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), scanTimeRange) / scanTimeRange);
        scores.push_back(pdScanTime); //6

        scores.push_back(std::max(candidateScores.cosineSim100MS1PreMono, 0.0)); //7
        scores.push_back(std::max(candidateScores.cosineSim100MS1Iso1, 0.0)); //7
        scores.push_back(std::max(candidateScores.cosineSim100MS1Iso2, 0.0)); //8
        scores.push_back(std::max(candidateScores.cosineSimSum45, 0.0)); //9
        scores.push_back(std::max(candidateScores.cosineSimSum20, 0.0)); //10
//        scores.push_back(std::max(candidateScores.cosineSim45MS1, 0.0)); //11
//        scores.push_back(std::max(candidateScores.cosineSim20MS1, 0.0)); //12

        const double pepLength = (-10.0 + candidateScores.peptideStringWithMods.size()) / 10.0;

        scores.push_back(pepLength); //13
        scores.push_back(candidateScores.theoFragmentCount); //14

        const QVector<double> cosineSimToAnchors
                = extractScoresFromVecFeatures(candidateScores.cosineSimToAnchorVec, theoMzIonsSize);
        scores.append(cosineSimToAnchors); //28-33

        const QVector<double> topHalfCosineSimScores = cosineSimToAnchors.mid(0, theoMzIonsSize / 2);
        const double topHalfCosineSimScoresSum = std::accumulate(topHalfCosineSimScores.begin(), topHalfCosineSimScores.end(), 0.0);
        scores.push_back(topHalfCosineSimScoresSum); //34

        const QVector<double> bottomHalfCosineSimScores = cosineSimToAnchors.mid(theoMzIonsSize / 2, theoMzIonsSize / 2);
        const double bottomHalfCosineSimScoresSum = std::accumulate(bottomHalfCosineSimScores.begin(), bottomHalfCosineSimScores.end(), 0.0);
        scores.push_back(bottomHalfCosineSimScoresSum); //35

        const double topBottomRatio
                = std::log(std::max(1.0, topHalfCosineSimScoresSum) / (topHalfCosineSimScoresSum + bottomHalfCosineSimScoresSum + 1.0));
        scores.push_back(topBottomRatio); //36

        const QVector<double> theoApexIntensity
                = extractScoresFromVecFeatures(candidateScores.theoIntensityVec, theoMzIonsSize);
        scores.append(theoApexIntensity);

        const QVector<double> intensityFoundMaxVec
                = extractScoresFromVecFeatures(candidateScores.intensityFoundMaxVec, theoMzIonsSize);
        const double totalIntensity = std::accumulate(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end(), 0.0001);
        const double totalIntensityLog = std::log(std::max(totalIntensity, std::numeric_limits<double>::min()));
        scores.push_back(totalIntensityLog);

        if (useNeuralNetworkScores) {

            const QVector<double> cosineSimShadowsToAnchorVec
                    = extractScoresFromVecFeatures(candidateScores.cosineSimShadowsToAnchorVec, theoMzIonsSize);
            scores.append(cosineSimShadowsToAnchorVec);

            scores.push_back(candidateScores.shadowsCosineSimSum);

            const QVector<double> shadowsIntensityRatioVec
                    = extractScoresFromVecFeatures(candidateScores.shadowsIntensityRatioVec, theoMzIonsSize);
            scores.append(shadowsIntensityRatioVec);

            for (double intzFound : intensityFoundMaxVec) {
                scores.push_back(intzFound / totalIntensity);
            }

            const int topSix = std::min(6, candidateScores.mzSearchedVec.size());
            QVector<double> ppmVec;
            for (int i = 0; i < topSix; i++) {

                const double mzSearched = candidateScores.mzSearchedVec.at(i);
                if (i >= candidateScores.mzFoundMeanVec.size()) {
                    break;
                }

                const double mzFound = candidateScores.mzFoundMeanVec[i];
                const double ppm = cosineSimToAnchors.at(i) * std::abs(mzFound - mzSearched) / mzSearched;
                ppmVec.push_back(ppm);
            }
        }

        const double charge = -2.0 + candidateScores.charge;
        scores.push_back(charge);
    }

    if (useNeuralNetworkScores) {
//        std::pow(std::max(0.0, candidateScores.cosineSimSpectrum), 3), //3
        scores.push_back(std::max(0.0, candidateScores.cosineSimSpectrum));
        scores.push_back(candidateScores.discriminateScore);
        scores.push_back(std::pow(std::max(0.0, candidateScores.klDivSpectrum), 3));

        QMap<QChar, int> aminoAcidCounts = {
                {'A', 0},
                {'C', 0},
                {'D', 0},
                {'E', 0},
                {'F', 0},
                {'G', 0},
                {'H', 0},
                {'I', 0},
                {'K', 0},
                {'L', 0},
                {'M', 0},
                {'N', 0},
                {'P', 0},
                {'Q', 0},
                {'R', 0},
                {'S', 0},
                {'T', 0},
                {'V', 0},
                {'W', 0},
                {'Y', 0}
        };

        for (const QChar aminoAcid : candidateScores.peptideStringWithMods) {

            if (!aminoAcidCounts.contains(aminoAcid)) {
                qDebug() << candidateScores.peptideStringWithMods << "missing amino acid";
            }

            e = ErrorUtils::isTrue(aminoAcidCounts.contains(aminoAcid)); ree;
            aminoAcidCounts[aminoAcid]++;
        }

        for (int cnt : aminoAcidCounts.values()) {
            scores.push_back(static_cast<double>(cnt));
        }

        const double mz = BiophysicalCalcs::calculateThomsonFromMass(
                candidateScores.mass,
                candidateScores.charge
        );
        scores.push_back((mz - 600.0) * 0.002);

        scores.push_back(candidateScores.peakShapeRatio1);
        scores.push_back(candidateScores.peakShapeRatio2);
        scores.push_back(candidateScores.peakShapeRatio3);
        scores.push_back(candidateScores.cosineSimSumBottom6);
        scores.push_back(candidateScores.cosineSimSumBottom6 / candidateScores.theoFragmentCount);

    }

    *scoreVec = scores;

    ERR_RETURN
}