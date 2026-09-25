#ifndef FRAGMENTCOMPETITIONTESTDATA_H
#define FRAGMENTCOMPETITIONTESTDATA_H

#include "CandidateScores.h"
#include "FragLibReaderRow.h"
#include "TargetDecoyCandidatePair.h"

// Owning fixtures must stay at a stable address because scores contain pointers.
struct CompetitionCandidateFixture {
    FragLibReaderRow library;
    TargetDecoyCandidatePair pair;
    CandidateScores scores;

    explicit CompetitionCandidateFixture(
        const QString &peptide = "PEPTIDEK",
        const QVector<float> &mzs = {300, 400, 500, 600, 700})
        : pair(PeptideStringWithMods(peptide), 0.0f) {
        library.mzVals = mzs;
        library.mass = 1000.0;
        library.precursorCharge = 2;
        QStringList labels;
        for (int ion = 0; ion < mzs.size(); ++ion) {
            library.intensityVals.push_back(1000.0f - ion);
            // Precursor-style test labels avoid introducing synthetic decoy
            // m/z shifts when a test isolates the target/decoy tie-breaker.
            labels.push_back("p");
        }
        library.ionLabels = labels.join(S_GLOBAL_SETTINGS.SEPARATOR);
        pair.setFragLibReaderRowPntr(&library);
        scores.targetDecoyCandidatePair = &pair;
        scores.featuresArray.fill(0.0f, FeaturesSize);
        scores.integrations.fill(0.0f, 12);
        scores.ionLabels.fill("p", 12);
        scores.featuresArray[Mass] = 1000.0f;
        scores.featuresArray[Charge] = 2.0f;
        scores.featuresArray[CosineSimSum100] = 6.0f;
        scores.scanNumber = 100;
        scores.scanTime = 60.0f;
        scores.scanTimeStart = 50.0f;
        scores.scanTimeEnd = 70.0f;
        scores.discriminantScore = 1.0;
        for (int ion = 0; ion < std::min(12, mzs.size()); ++ion) {
            scores.integrations[ion] = 100.0f;
            scores.featuresArray[CosineSimToAnchor1 + ion] = 1.0f;
            scores.featuresArray[MzFoundMean1 + ion] = mzs.at(ion);
            scores.featuresArray[IntensityFoundMax1 + ion] = 1.0f;
        }
    }

    CompetitionCandidateFixture(const CompetitionCandidateFixture &) = delete;
    CompetitionCandidateFixture &operator=(const CompetitionCandidateFixture &) = delete;
};

#endif
