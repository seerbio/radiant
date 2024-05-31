
#include "MsReaderMzMLMapped.h"

#include "ErrorUtils.h"
#include "ParallelUtils.h"
#include "SqlUtils.h"

#include <QCoreApplication>
#include <QFile>
#include <QMap>
#include <QXmlStreamReader>

#include <zlib.h>

#include <algorithm>
#include <iostream>
#include <QDebug>
#include <QFileInfo>


const QString scanKey = QStringLiteral("scan");

const QString dataProcessingElementName = QStringLiteral("<dataProcessing");
const QString mzmlElementName = QStringLiteral("<mzML ");
const QString offsetElementName = QStringLiteral("<offset ");
const QString peaksElementName = QStringLiteral("<peaks ");
const QString spectrumElementName = QStringLiteral("<spectrum ");
const QString spectrumElementEndName = QStringLiteral("</spectrum");
const QString binaryElementName = QStringLiteral("<binary>");
const QString binaryElementEndName = QStringLiteral("</binary>");

const QString BASEPEAK_INTENSITY = QStringLiteral("MS:1000505");
const QString BINARY_DATA_ARRAY = QStringLiteral("binaryDataArray");
const QString BINARY_DATA_ARRAY_LIST = QStringLiteral("binaryDataArrayList");
const QString BINARY = QStringLiteral("binary");
const QString CHARGE_STATE = QStringLiteral("MS:1000041");
const QString COLLISION_ENERGY = QStringLiteral("MS:1000045");
const QString ACCESSION = QStringLiteral("accession");
const QString NAME = QStringLiteral("name");
const QString VALUE = QStringLiteral("value");
const QString BIT64_FLOAT = QStringLiteral("MS:1000523");
const QString BIT32_FLOAT = QStringLiteral("MS:1000521");
const QString FILTER_STRING = QStringLiteral("MS:1000512");
const QString FAIMS_VOLTAGE = QStringLiteral("MS:1001581");
const QString INTENSITY_ARRAY = QStringLiteral("MS:1000515");
const QString NO_COMPRESSION = QStringLiteral("MS:1000576");
const QString MS_LEVEL = QStringLiteral("MS:1000511");
const QString MZ_ARRAY = QStringLiteral("MS:1000514");
const QString PRECURSOR = QStringLiteral("precursor");
const QString SCAN = QStringLiteral("scan");
const QString SCAN_START_TIME = QStringLiteral("MS:1000016");
const QString TIC = QStringLiteral("MS:1000285");
const QString UNIT_ACCESSION = QStringLiteral("unitAccession");
const QString ZLIB_COMPRESSION = QStringLiteral("MS:1000574");

const QString MILLISECOND = QStringLiteral("UO:0000028");
const QString SECOND = QStringLiteral("UO:0000010");
const QString MINUTE = QStringLiteral("UO:0000031");
const QString MZ = QStringLiteral("MS:1000040");

const QString  PRECURSOR_TARGET_MZ = QStringLiteral("MS:1000827");
const QString  PRECURSOR_LOWER_WINDOW_OFFSET = QStringLiteral("MS:1000828");
const QString  PRECURSOR_UPPER_WINDOW_OFFSET = QStringLiteral("MS:1000829");


enum TYPES {
    FLOAT32 = 0,
    FLOAT64,
    NULL_TYPE
};


enum COMPRESSION {
    ZLIB = 0,
    NO_COMPRESS,
    NULL_COMPRESSION
};


enum SPECTRUM_TYPE {
    MZ_ARR = 0,
    INTENSITY_ARR,
    NULL_SPECTRUM_TYPE
};


class Q_DECL_HIDDEN MsReaderMzMLMapped::PrivateData  {

public:

    PrivateData(
            QMap<ScanNumber, MsScanInfo> *msScanInfo,
            QMap<ScanNumber, ScanPoints> *scanPoints
                );

    ~PrivateData();

    Err openFile(const QString& filename);
    Err closeFile();

    [[nodiscard]] bool isScanNumberValid(int scanNumber) const;

public:
    QMap<ScanNumber, MsScanInfo> *m_msScanInfo;
    QMap<ScanNumber, ScanPoints> *m_scanPoints;

};

MsReaderMzMLMapped::PrivateData::PrivateData(
        QMap<ScanNumber, MsScanInfo> *msScanInfo,
        QMap<ScanNumber, ScanPoints> *scanPoints
        )
: m_msScanInfo(msScanInfo)
, m_scanPoints(scanPoints) {}

MsReaderMzMLMapped::PrivateData::~PrivateData() {
    closeFile();
}

namespace {

    struct FileChunk {
        const char* data;
        size_t start;
        size_t end;
        size_t overlap;
    };

    QVector<FileChunk> buildFileChunks(
        qint64 fileSize,
        uchar* ucharData
        ) {

        QVector<FileChunk> fileChunks;

        const auto* data = reinterpret_cast<const char*>(ucharData);
        const int numChunks = ParallelUtils::numberOfAvailableSystemProcessors();

        const size_t overlapSize = 1024 * 500; // 500 KB overlap
        const size_t chunkSize = fileSize / numChunks;

        for (int i = 0; i < numChunks; ++i) {
            size_t start = i * chunkSize;
            size_t end = (i == numChunks - 1) ? fileSize : (i + 1) * chunkSize + overlapSize;
            end = std::min(end, static_cast<size_t>(fileSize));

            if (i != 0) {
                while (start > 0 && data[start] != '>' && data[start - 1] != '<') {
                    --start;
                }
            }
            if (i != numChunks - 1) {
                while (end < fileSize && data[end] != '>' && data[end - 1] != '<') {
                    ++end;
                }
            }

            fileChunks.append({data, start, end, overlapSize});
        }

        return fileChunks;
    }

    bool isCurrentTag(const QXmlStreamReader &reader, const QString &tagName) {
        return (0 == QStringRef::compare(reader.name(), tagName, Qt::CaseInsensitive));
    }

    bool decompress(const QByteArray &input, QByteArray *output) {

        Q_ASSERT(output);

        const size_t gzipWindowsBits = 15;
        const uInt gzipChunkSize = 32 * 1024;

        output->clear();

        if (input.isEmpty()) {
            return false;
        }

        z_stream strm;
        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;
        strm.avail_in = 0;
        strm.next_in = Z_NULL;

        if (Z_OK != inflateInit2(&strm, gzipWindowsBits)) {
            return false;
        }

        const char *inputData = input.data();
        uInt inputDataLeft = input.length();

        int ret = Z_OK;

        do {

            const uInt chunkSize = std::min(gzipChunkSize, inputDataLeft);

            if (chunkSize <= 0) {
                break;
            }

            strm.next_in = reinterpret_cast<z_const Bytef *>(const_cast<char *>(inputData));
            strm.avail_in = chunkSize;

            inputData += chunkSize;
            inputDataLeft -= chunkSize;

            do {

                char out[gzipChunkSize];

                strm.next_out = reinterpret_cast<Bytef *>(out);
                strm.avail_out = gzipChunkSize;

                ret = inflate(&strm, Z_NO_FLUSH);

                switch (ret) {
                    case Z_NEED_DICT:
                        ret = Z_DATA_ERROR;
                    case Z_DATA_ERROR:
                    case Z_MEM_ERROR:
                    case Z_STREAM_ERROR:
                        inflateEnd(&strm);
                        return false;
                }

                const int have = (gzipChunkSize - strm.avail_out);

                if (have > 0) {
                    output->append(out, have);
                }
            } while (strm.avail_out == 0);
        } while (ret != Z_STREAM_END);

        inflateEnd(&strm);

        return ret == Z_STREAM_END;
    }

    Err getScanNumber(const QXmlStreamReader &reader, int *scanNumber) {

        ERR_INIT

        e = ErrorUtils::isTrue(isCurrentTag(reader, spectrumElementName)); ree;

        const QXmlStreamAttributes attributesSpectrumTag = reader.attributes();

        bool ok = true;
        const QString id = attributesSpectrumTag.value("id").toString();

        *scanNumber = id.split("=", Qt::SkipEmptyParts).back().toInt(&ok);
        e = ErrorUtils::isTrue(ok); ree;

        ERR_RETURN
    }

    bool stringMatch(const QString &s1, const QString &s2) {
        return (0 == QString::compare(s1, s2, Qt::CaseInsensitive) );
    }

    bool stringMatch(const QStringRef &s1, const QString &s2) {
        return (0 == QStringRef::compare(s1, s2, Qt::CaseInsensitive) );
    }

    void filterZeroIntensityQPoints(ScanPoints *vec) {

        const auto terminatorLogic = [](const ScanPoint &p){
            return MathUtils::tZero(p.y());
        };

        const auto terminator = std::remove_if(vec->begin(), vec->end(), terminatorLogic);
        vec->erase(terminator, vec->end());
    }

    Err processBinaryData(QXmlStreamReader &reader, ScanPoints *scanPoints) {

        ERR_INIT

        TYPES type = NULL_TYPE;
        COMPRESSION compression = COMPRESSION::NULL_COMPRESSION;
        SPECTRUM_TYPE spectrumType = SPECTRUM_TYPE::NULL_SPECTRUM_TYPE;

        QVector<double> mzValues;
        QVector<double> intensityValues;

        while (!reader.atEnd() && !reader.hasError()) {

            reader.readNext();

            if (reader.isStartElement()) {

                const QXmlStreamAttributes attributesCurrentTag = reader.attributes();
                const QStringRef currentAccession = attributesCurrentTag.value(ACCESSION);
                const QStringRef currentReaderName = reader.name();

                if (stringMatch(currentAccession, BIT64_FLOAT)) {
                    type = TYPES::FLOAT64;
                }
                if (stringMatch(currentAccession, BIT32_FLOAT)) {
                    type = TYPES::FLOAT32;
                }
                if (stringMatch(currentAccession, ZLIB_COMPRESSION)) {
                    compression = COMPRESSION::ZLIB;
                }
                if (stringMatch(currentAccession, NO_COMPRESSION)) {
                    compression = COMPRESSION::NO_COMPRESS;
                }
                if (stringMatch(currentAccession, MZ_ARRAY)) {
                    spectrumType = SPECTRUM_TYPE::MZ_ARR;
                }
                if (stringMatch(currentAccession, INTENSITY_ARRAY)) {
                    spectrumType = SPECTRUM_TYPE::INTENSITY_ARR;
                }

                if (stringMatch(currentReaderName, BINARY)) {

                    const std::string elementText = reader.readElementText().toStdString();
                    QByteArray binaryData = QByteArray::fromBase64(QByteArray::fromStdString(elementText));

                    QByteArray binaryDataOut;
                    if (compression == COMPRESSION::ZLIB) {
                        decompress(binaryData, &binaryDataOut);
                    }
                    else if (compression == COMPRESSION::NO_COMPRESS){
                        binaryDataOut = binaryData;
                    }

                    QVector<double> decodedBLOB;
                    if (type == TYPES::FLOAT64) {
                        decodedBLOB = SqlUtils::decodeBLOB<double>(binaryDataOut);
                    }
                    else if (type == TYPES::FLOAT32) {
                        const QVector<float> decodedBLOBFloat = SqlUtils::decodeBLOB<float>(binaryDataOut);
                        decodedBLOB.reserve(decodedBLOBFloat.size());
                        std::copy(decodedBLOBFloat.cbegin(), decodedBLOBFloat.cend(), std::back_inserter(decodedBLOB));
                    }

                    if (spectrumType == SPECTRUM_TYPE::MZ_ARR) {
                        mzValues = decodedBLOB;
                    }
                    else if (spectrumType == SPECTRUM_TYPE::INTENSITY_ARR) {
                        intensityValues = decodedBLOB;
                    }
                }

                continue;
            }

            if (reader.isEndElement()) {

                if (isCurrentTag(reader, BINARY_DATA_ARRAY_LIST)) {

                    e = ErrorUtils::isEqual(mzValues.size(), intensityValues.size()); ree;

                    scanPoints->resize(mzValues.size());
                    scanPoints->reserve(mzValues.size());

                    for (int i = 0; i < mzValues.size(); i++) {
                        const ScanPoint point(
                                static_cast<float>(mzValues.at(i)),
                                static_cast<float>(intensityValues.at(i))
                                );
                        (*scanPoints)[i] = point;
                    }

                    filterZeroIntensityQPoints(scanPoints);

                    ERR_RETURN
                }
            }
        }

        ERR_RETURN
    }

    Err processScanData(QXmlStreamReader &reader, MsScanInfo *msScanInfo) {

        ERR_INIT

        while (!reader.atEnd() && !reader.hasError()) {

            reader.readNext();

            if (reader.isStartElement()) {

                const QXmlStreamAttributes attributesCurrentTag = reader.attributes();
                const QStringRef currentAccession = attributesCurrentTag.value(ACCESSION);

                if (stringMatch(currentAccession, SCAN_START_TIME)) {

                    bool ok = true;
                    msScanInfo->scanTime = attributesCurrentTag.value(VALUE).toDouble(&ok);
                    e = ErrorUtils::isTrue(ok); ree;

                }

                continue;
            }

            if (reader.isEndElement()) {
                if (isCurrentTag(reader, SCAN)) {
                    ERR_RETURN
                }
            }
        }

        ERR_RETURN
    }

    Err processPrecursor(QXmlStreamReader &reader, MsScanInfo *msScanInfo) {

        ERR_INIT

        while (!reader.atEnd() && !reader.hasError()) {

            reader.readNext();

            if (reader.isStartElement()) {

                const QXmlStreamAttributes attributesCurrentTag = reader.attributes();
                const QStringRef currentAccession = attributesCurrentTag.value(ACCESSION);

                if (stringMatch(currentAccession, PRECURSOR_TARGET_MZ)) {
                    bool ok = true;
                    msScanInfo->precursorTargetMz = attributesCurrentTag.value(VALUE).toDouble(&ok);
                    e = ErrorUtils::isTrue(ok); ree;
                }

                if (stringMatch(currentAccession, PRECURSOR_LOWER_WINDOW_OFFSET)) {
                    bool ok = true;
                    msScanInfo->isoWindowLower = attributesCurrentTag.value(VALUE).toDouble(&ok);
                    e = ErrorUtils::isTrue(ok); ree;
                }

                if (stringMatch(currentAccession, PRECURSOR_UPPER_WINDOW_OFFSET)) {
                    bool ok = true;
                    msScanInfo->isoWindowUpper = attributesCurrentTag.value(VALUE).toDouble(&ok);
                    e = ErrorUtils::isTrue(ok); ree;
                }

                if (stringMatch(currentAccession, COLLISION_ENERGY)) {
                    bool ok = true;
                    msScanInfo->collisionEnergy = static_cast<int>(attributesCurrentTag.value(VALUE).toDouble(&ok));
                    e = ErrorUtils::isTrue(ok); ree;
                }

                continue;
            }

            if (reader.isEndElement()) {
                if (isCurrentTag(reader, PRECURSOR)) {
                    ERR_RETURN
                }
            }
        }

        ERR_RETURN
    }

    QPair<Err, QStringList> parseScan(QXmlStreamReader &reader) {

        ERR_INIT

        MsScanInfo msScanInfo;

        QStringList returnList;

        e = getScanNumber(reader, &msScanInfo.scanNumber); rree;

        while (!reader.atEnd() && !reader.hasError()) {

            reader.readNext();

            if (reader.isStartElement()) {

                const QXmlStreamAttributes attributesCurrentTag = reader.attributes();
                const QStringRef currentAccession = attributesCurrentTag.value(ACCESSION);
                const QStringRef currentValue = attributesCurrentTag.value(VALUE);
                const QStringRef currentReaderName = reader.name();

                bool ok = true;

                if (stringMatch(currentAccession, MS_LEVEL)) {
                    msScanInfo.msLevel = currentValue.toInt(&ok);
                    e = ErrorUtils::isTrue(ok); rree;
                }

    //            if (stringMatch(currentAccession, FAIMS_VOLTAGE)) {
    //                msScanInfo.faimsVoltage = static_cast<int>(currentValue.toDouble(&ok));
    //                e = ErrorUtils::isTrue(ok); ree;
    //            }

    //            if (stringMatch(currentAccession, TIC)) {
    //                msScanInfo.TIC = currentValue.toDouble(&ok);
    //                e = ErrorUtils::isTrue(ok); ree;
    //            }
    //
    //            if (stringMatch(currentAccession, BASEPEAK_INTENSITY)) {
    //                msScanInfo.basePeakIntensity = currentValue.toDouble(&ok);
    //                e = ErrorUtils::isTrue(ok); ree;
    //            }

                if (stringMatch(currentReaderName, SCAN)) {
                    e = processScanData(reader, &msScanInfo);  rree;
                }

                if (stringMatch(currentReaderName, PRECURSOR)) {
                    e = processPrecursor(reader, &msScanInfo); rree;
                }

                if (stringMatch(currentReaderName, BINARY_DATA_ARRAY)) {
                    QVector<ScanPoint> scanPoints;
                    e = processBinaryData(reader, &scanPoints); rree;
                    // (*m_scanPoints)[msScanInfo.scanNumber] = scanPoints;
                }

                continue;
            }

            if (reader.isEndElement()) {
                if (isCurrentTag(reader, spectrumElementName)) {

                    // m_msScanInfo->insert(msScanInfo.scanNumber, msScanInfo);

                    returnList.push_back(QString::number(msScanInfo.scanNumber));

                    if (msScanInfo.scanNumber % 1000 == 0) {
                        qDebug() << "mzML Scans Parsed" << msScanInfo.scanNumber;
                    }

                    return {e, returnList};
                }
            }
        }

        return {eError, {}};
    }

    QPair<Err, QStringList> parseMsRun(QXmlStreamReader &reader) {

        ERR_INIT

        const QString SCAN_COUNT_ATTRIBUTE_NAME = QStringLiteral("count");

        bool ok = true;

        QXmlStreamAttributes attributes = reader.attributes();
        const int scanCount = attributes.value(SCAN_COUNT_ATTRIBUTE_NAME).toInt(&ok);
        e = ErrorUtils::isTrue(ok); rree;

        qDebug() << "Scan Count:" << scanCount;

        QStringList returnLister;

        while (!reader.atEnd()) {

            reader.readNext();

            if (reader.isStartElement()) {
                if (isCurrentTag(reader, spectrumElementName)) {
                    const QPair<Err, QStringList> parsedScan = parseScan(reader);
                    e = parsedScan.first;
                    for (const QString s : parsedScan.second) {
                        returnLister.push_back(s);
                    }
                    continue;
                }
            }
            if (reader.isEndElement()) {
                if (isCurrentTag(reader, mzmlElementName)) {
                    break;
                }
            }
        }

        // e = ErrorUtils::isEqual(scanCount, m_msScanInfo->size()); eee_absorb;
        return {e, returnLister};
    }

    QMap<QString, QString> parseAttributes(const QString& xmlString) {
        QMap<QString, QString> attributesMap;

        QXmlStreamReader xml(xmlString);
        bool elementProcessed = false;

        while (!xml.atEnd() && !xml.hasError()) {
            QXmlStreamReader::TokenType token = xml.readNext();
            if (token == QXmlStreamReader::StartElement) {
                elementProcessed = true;
                foreach (const QXmlStreamAttribute &attr, xml.attributes()) {
                    attributesMap.insert(attr.name().toString(), attr.value().toString());
                }
            }
        }

        // Suppress the premature end of document warning if we successfully processed the element
        if (xml.hasError() && !(xml.error() == QXmlStreamReader::PrematureEndOfDocumentError && elementProcessed)) {
            qWarning() << "XML error:" << xml.errorString();
        }

        return attributesMap;
    }

    Err parseStringForScanNumber(
        const QString& str,
        ScanNumber *scanNumber
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(str); ree;

        QHash<QString, int> result;
        QStringList pairs = str.split(' ');

        for(const QString& pair : pairs) {
            QStringList keyValue = pair.split('=');
            if(keyValue.size() == 2) {
                bool ok;
                result[keyValue[0]] = keyValue[1].toInt(&ok);
                e = ErrorUtils::isTrue(ok); ree;
            }
        }

        *scanNumber = result.value(scanKey);
        ERR_RETURN
    }

    Err processSpectrumKey(
        const QString &str,
        ScanNumber *scanNumber
        ) {

        ERR_INIT

        const QMap<QString, QString> m = parseAttributes(str.trimmed());

        try {
            e = parseStringForScanNumber(
                m.value("id"),
                scanNumber
                );

            try {
                e = ErrorUtils::toInt(m.value("index"), scanNumber);
                *scanNumber += 1;
            }
            catch (std::exception &ex) {
                qCritical() << QString::fromStdString(ex.what());
                rrr(eValueError);
            }
        }
        catch (std::exception &ex) {
            qWarning() << QString::fromStdString(ex.what());
            rrr(eValueError);
        }

        ERR_RETURN
    }

    QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>> processChunk(const FileChunk& chunk) {

        ERR_INIT

        QVector<QPair<MsScanInfo, ScanPoints>> msScanInfosVsScanPoints;

        MsScanInfo msScanInfoLocal;
        ScanPoints scanPointsLocal;

        TYPES type = NULL_TYPE;
        COMPRESSION compression = COMPRESSION::NULL_COMPRESSION;
        SPECTRUM_TYPE spectrumType = SPECTRUM_TYPE::NULL_SPECTRUM_TYPE;

        QVector<double> mzValues;
        QVector<double> intensityValues;

        QString str;
        for (size_t i = chunk.start; i < chunk.end; ++i) {

            str += chunk.data[i];

            if (chunk.data[i] == '\n') {

                str = str.trimmed();
                if (str.isEmpty() || str.front() != '<') {
                    str.clear();
                    continue;
                }

                if (str.contains(spectrumElementName)) {
                    e = processSpectrumKey(str, &msScanInfoLocal.scanNumber); rree;
                }
                else if (str.contains(MS_LEVEL)) {
                    QMap<QString, QString> attributes = parseAttributes(str);
                    e = ErrorUtils::toInt(
                        attributes.value("value"),
                        &msScanInfoLocal.msLevel
                    ); rree;
                }
                else if (str.contains(SCAN_START_TIME)) {
                    QMap<QString, QString> attributes = parseAttributes(str);

                    e = ErrorUtils::toFloat(
                        attributes.value("value"),
                        &msScanInfoLocal.scanTime
                    ); rree

                    msScanInfoLocal.scanTime = attributes.value("unitName").contains("minute")
                                             ? msScanInfoLocal.scanTime
                                             : msScanInfoLocal.scanTime / 60.0f;
                }
                else if (str.contains(PRECURSOR_TARGET_MZ)) {
                    QMap<QString, QString> attributes = parseAttributes(str);
                    e = ErrorUtils::toFloat(
                        attributes.value("value"),
                        &msScanInfoLocal.precursorTargetMz
                    ); rree;
                }
                else if (str.contains(PRECURSOR_LOWER_WINDOW_OFFSET)) {
                    QMap<QString, QString> attributes = parseAttributes(str);
                    e = ErrorUtils::toFloat(
                        attributes.value("value"),
                        &msScanInfoLocal.isoWindowLower
                    ); rree;
                }
                else if (str.contains(PRECURSOR_UPPER_WINDOW_OFFSET)) {
                    QMap<QString, QString> attributes = parseAttributes(str);
                    e = ErrorUtils::toFloat(
                        attributes.value("value"),
                        &msScanInfoLocal.isoWindowUpper
                    ); rree;
                }
                else if (str.contains(COLLISION_ENERGY)) {
                    QMap<QString, QString> attributes = parseAttributes(str);
                    e = ErrorUtils::toFloat(
                        attributes.value("value"),
                        &msScanInfoLocal.collisionEnergy
                    ); rree;
                }
                else if (str.contains(BIT64_FLOAT)) {
                    type = TYPES::FLOAT64;
                }
                else if (str.contains(BIT32_FLOAT)) {
                    type = TYPES::FLOAT32;
                }
                else if (str.contains(ZLIB_COMPRESSION)) {
                    compression = COMPRESSION::ZLIB;
                }
                else if (str.contains(NO_COMPRESSION)) {
                    compression = COMPRESSION::NO_COMPRESS;
                }
                else if (str.contains(MZ_ARRAY)) {
                    spectrumType = SPECTRUM_TYPE::MZ_ARR;
                }
                else if (str.contains(INTENSITY_ARRAY)) {
                    spectrumType = SPECTRUM_TYPE::INTENSITY_ARR;
                }
                else if (str.contains(binaryElementName)) {

                    str = str.replace(binaryElementName, "");
                    str = str.replace(binaryElementEndName, "");
                    const std::string elementText = str.toStdString();

                    QByteArray binaryData = QByteArray::fromBase64(QByteArray::fromStdString(elementText));

                    QByteArray binaryDataOut;
                    if (compression == COMPRESSION::ZLIB) {
                        decompress(binaryData, &binaryDataOut);
                    }
                    else if (compression == COMPRESSION::NO_COMPRESS){
                        binaryDataOut = binaryData;
                    }

                    QVector<double> decodedBLOB;
                    if (type == TYPES::FLOAT64) {
                        decodedBLOB = SqlUtils::decodeBLOB<double>(binaryDataOut);
                    }
                    else if (type == TYPES::FLOAT32) {
                        const QVector<float> decodedBLOBFloat = SqlUtils::decodeBLOB<float>(binaryDataOut);
                        decodedBLOB.reserve(decodedBLOBFloat.size());
                        std::copy(decodedBLOBFloat.cbegin(), decodedBLOBFloat.cend(), std::back_inserter(decodedBLOB));
                    }

                    if (spectrumType == SPECTRUM_TYPE::MZ_ARR) {
                        mzValues = decodedBLOB;
                    }
                    else if (spectrumType == SPECTRUM_TYPE::INTENSITY_ARR) {
                        intensityValues = decodedBLOB;
                    }
                }
                else if (str.contains(spectrumElementEndName)) {

                    if (mzValues.size() != intensityValues.size() || msScanInfoLocal.scanNumber < 0) {
                        msScanInfoLocal = MsScanInfo();
                        scanPointsLocal = ScanPoints();
                        str.clear();
                        continue;
                    }

                    e = ErrorUtils::isEqual(mzValues.size(), intensityValues.size()); rree;

                    scanPointsLocal.resize(mzValues.size());
                    scanPointsLocal.reserve(mzValues.size());

                    for (int i = 0; i < mzValues.size(); i++) {
                        const ScanPoint point(
                                static_cast<float>(mzValues.at(i)),
                                static_cast<float>(intensityValues.at(i))
                                );
                        scanPointsLocal[i] = point;
                    }

                    msScanInfosVsScanPoints.push_back({msScanInfoLocal, scanPointsLocal});

                    msScanInfoLocal = MsScanInfo();
                    scanPointsLocal = ScanPoints();
                }

                str.clear();
            }
        }

        return {e, msScanInfosVsScanPoints};
    }

}//NAMESPACE
Err MsReaderMzMLMapped::PrivateData::openFile(const QString &filename) {

    ERR_INIT

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        qCritical() << "Error opening file:" << file.errorString();
        ERR_RETURN
    }

    const QFileInfo fileInfo(file);
    const qint64 fileSize = fileInfo.size();
    uchar* ucharData = file.map(0, fileSize);
    if (!ucharData) {
        qCritical() << "Error mapping file:" << file.errorString();
        ERR_RETURN
    }

    const QVector<FileChunk> chunks = buildFileChunks(fileSize, ucharData);

#define RUN_PARALLEL
#ifdef RUN_PARALLEL
    QFuture<QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>>> future = QtConcurrent::mapped(chunks, processChunk);
    future.waitForFinished();

    for (const QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>> &result : future) {
        e = result.first; ree;
        for (const QPair<MsScanInfo, ScanPoints> &pr : result.second) {
            const MsScanInfo &msScanInfo = pr.first;
            const ScanPoints &scanPoints = pr.second;
            if (m_msScanInfo->contains(msScanInfo.scanNumber) || msScanInfo.scanNumber < 0) {
                continue;
            }

            m_msScanInfo->insert(msScanInfo.scanNumber, msScanInfo);
            m_scanPoints->insert(msScanInfo.scanNumber, scanPoints);
        }
    }
    qDebug() << "Scan Count" << m_msScanInfo->size();
#else

    for (const FileChunk &fc : chunks) {
        QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>> result = processChunk(fc);
        e = result.first; ree;
        for (const QPair<MsScanInfo, ScanPoints> &pr : result.second) {
            const MsScanInfo &msScanInfo = pr.first;
            const ScanPoints &scanPoints = pr.second;
            if (m_msScanInfo->contains(msScanInfo.scanNumber) || msScanInfo.scanNumber < 0) {
                continue;
            }
            m_msScanInfo->insert(msScanInfo.scanNumber, msScanInfo);
            m_scanPoints->insert(msScanInfo.scanNumber, scanPoints);
        }
    }

#endif

    file.unmap(ucharData);
    file.close();

    ERR_RETURN
}

Err MsReaderMzMLMapped::PrivateData::closeFile() {

    ERR_INIT

    m_msScanInfo->clear();
    m_scanPoints->clear();

    ERR_RETURN
}

bool MsReaderMzMLMapped::PrivateData::isScanNumberValid(int scanNumber) const {
    return scanNumber >= 1 && scanNumber <= m_msScanInfo->size();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

MsReaderMzMLMapped::MsReaderMzMLMapped() {
    m_d.reset(new PrivateData(
            &this->m_msScanInfo,
            &this->m_scanPoints)
            );
}

MsReaderMzMLMapped::~MsReaderMzMLMapped() {
}

Err MsReaderMzMLMapped::openFile(const QString &filePath) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(filePath); ree;
    m_filePath = filePath;

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            MzMLNamespace::MZML_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

    e = m_d->openFile(filePath); ree;

    e = printFileInfo(); ree;
    ERR_RETURN
}

Err MsReaderMzMLMapped::closeFile() {
    return m_d->closeFile();
}

MsReaderBase MsReaderMzMLMapped::msReaderBase() {

    MsReaderBase msReaderBase;
    // NOTE: it is important to set scanInfos first or setScanPoints will Err out.
    msReaderBase.setMsScanInfo(*m_d->m_msScanInfo);
    msReaderBase.setScanPoints(*m_d->m_scanPoints);

    return msReaderBase;
}
