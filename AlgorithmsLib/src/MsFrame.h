 //
// Created by anichols on 4/8/23.
//

#ifndef PYTHIADIACPP_MSFRAME_H
#define PYTHIADIACPP_MSFRAME_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"
#include "PythiaParameterReader.h"

using namespace Error;

namespace MsFrameScanPointRowsNamespace {
    const QString FRAME_INDEX = QStringLiteral("frameIndex");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");

    const QStringList keysToCheck = {
         FRAME_INDEX,
         MZ_VALS,
         INTENSITY_VALS
    };
}//namespace
struct MsFrameScanPointRows : public ParquetReaderInputBase {

    FrameIndex frameIndex = -1;
    QVector<double> mzVals;
    QVector<double> intensityVals;

    QMap<QString, QVariant> map() override {

        using namespace MsFrameScanPointRowsNamespace;

        return {
            {FRAME_INDEX, QVariant(frameIndex)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MsFrameScanPointRowsNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frameIndex = dataMap.value(FRAME_INDEX).toInt();
        mzVals = bytesArrayToQVector<double>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());

        ERR_RETURN
    }

};


class ALGORITHMSLIB_EXPORTS MsFrame {

public:

    friend class MsFrameTests;

    MsFrame() = default;
    ~MsFrame() = default;

    Err init(
            const QMap<ScanNumber, ScanPoints> &scanPoints,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
            );

    [[nodiscard]] bool isValid() const ;

    static Err writeFrameScans(
            const QMap<FrameIndex, ScanPoints> &framesVsScanPoints,
            const QString &outputFilePath
            );

    [[nodiscard]] int scanCount() const;

    [[nodiscard]] QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints() const;

    [[nodiscard]] QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints() const;

    [[nodiscard]] ScanNumber scanNumberFromFrameIndex(FrameIndex frameIndex) const;

    [[nodiscard]] ScanTime scanTimeFromScanNumber(ScanNumber scanNumber) const;

    [[nodiscard]] ScanNumber scanNumberFromScanTime(ScanTime scanTime) const;

    [[nodiscard]] ScanNumber frameIndexFromScanNumber(ScanNumber scanNumber) const;

    [[nodiscard]] ScanPoints getScanPointsByScanNumber(ScanNumber scanNumber) const;

private:

    Err buildFrameIndexVsScanNumber();

private:

    QMap<ScanNumber, ScanPoints> m_frame;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QMap<FrameIndex, ScanNumber> m_frameIndexVsScanNumber;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

};


#endif //PYTHIADIACPP_MSFRAME_H
