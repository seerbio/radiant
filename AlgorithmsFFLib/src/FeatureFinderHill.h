//
// Created by anichols on 12/11/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILL_H
#define PYTHIACPP_FEATUREFINDERHILL_H

#include "AlgorithmsFFLib_Exports.h"
#include "BiophysicalCalcs.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

#include <QVector>


using namespace Error;

namespace FeatureFinderHillEntryNamespace {

    const QString SCAN_NUMBER_VALS = QStringLiteral("scanNumberVals");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString FRAME_INDEX_VALS = QStringLiteral("frameIndexVals");

    const QStringList keysToCheck = {
            SCAN_NUMBER_VALS,
            MZ_VALS,
            INTENSITY_VALS,
            FRAME_INDEX_VALS
    };

}

struct FILEREADERSLIB_EXPORTS FeatureFinderHillEntry : public ParquetReaderInputBase {

    QVector<ScanNumber> scanNumberVals;
    QVector<double> mzVals;
    QVector<float> intensityVals;
    QVector<FrameIndex> frameIndexVals;

    QMap<QString, QVariant> map() override {

        using namespace FeatureFinderHillEntryNamespace;

        return {
                {SCAN_NUMBER_VALS, QVariant(qVectorToQByteArray(scanNumberVals))},
                {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
                {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
                {FRAME_INDEX_VALS, QVariant(qVectorToQByteArray(frameIndexVals))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace FeatureFinderHillEntryNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        scanNumberVals = bytesArrayToQVector<ScanNumber>(dataMap.value(SCAN_NUMBER_VALS).toByteArray());
        mzVals = bytesArrayToQVector<double>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS).toByteArray());
        frameIndexVals = bytesArrayToQVector<FrameIndex>(dataMap.value(FRAME_INDEX_VALS).toByteArray());

        ERR_RETURN
    }

};


class ALGORITHMSFFLIB_EXPORTS FeatureFinderHill {

public:

    FeatureFinderHill() = default;
    ~FeatureFinderHill() = default;

    bool operator == (const FeatureFinderHill &c) const {
        if (mzMean() == c.mzMean() && maxIntensityScanNumberIndex() == c.maxIntensityScanNumberIndex())
            return true;
        return false;
    }

    void addPoint(
        ScanNumberIndex scanNumberIndex,
        ScanNumber scanNumber,
        double mzVal,
        float intensityVal
    );

    FeatureFinderHillEntry featureFinderHillEntry();

    Err initFromFeatureFinderHillEntry(const FeatureFinderHillEntry &featureFinderHillEntry);

    [[nodiscard]] double mzMean() const;

    [[nodiscard]] double mzStDev() const;

    [[nodiscard]] QPair<double, double> mzMinMax() const;

    [[nodiscard]] int maxIntensityScanNumber() const;

    [[nodiscard]] int maxIntensityScanNumberIndex() const;

    [[nodiscard]] float intensityValueMax() const;

    [[nodiscard]] int scanCount() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberMinMax() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberIndexMinMax() const;

    [[nodiscard]] QVector<int> scanNumbers() const;

    [[nodiscard]] QVector<int> scanNumberIndexes() const;

    [[nodiscard]] QVector<double> mzVals() const;

    [[nodiscard]] QVector<float> intensities() const;


private:

    QVector<double> m_mzVals;
    QVector<ScanNumberIndex> m_scanNumberIndexes;
    QVector<ScanNumber> m_scanNumbers;
    QVector<float> m_intensities;

};


#endif //PYTHIACPP_FEATUREFINDERHILL_H
