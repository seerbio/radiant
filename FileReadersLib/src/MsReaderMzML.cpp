
#include "MsReaderMzML.h"

#include "ErrorUtils.h"
#include "SqlUtils.h"

#include <QCoreApplication>
#include <QFile>
#include <QMap>
#include <QXmlStreamReader>

//Compiled using visual studio: /PythiaCpp/ThirdPartyLibs/src/zlib-1.2.11/contrib/vstudio/vc14/zlibvc.sln
#include <zlib.h>

#include <algorithm>
#include <iostream>
#include <QDebug>
#include <QFileInfo>


const QString dataProcessingElementName = QStringLiteral("dataProcessing");
const QString mzmlElementName = QStringLiteral("mzML");
const QString offsetElementName = QStringLiteral("offset");
const QString peaksElementName = QStringLiteral("peaks");
const QString spectrumElementName = QStringLiteral("spectrum");
const QString spectrumListElementName = QStringLiteral("spectrumList");

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


class Q_DECL_HIDDEN MsReaderMzML::PrivateData  {

public:

    PrivateData(
            QMap<ScanNumber, MsScanInfo> *msScanInfo,
            QMap<ScanNumber, ScanPoints> *scanPoints
                );

    ~PrivateData();

    Err openFile(const QString& filename);
    Err closeFile();

    [[nodiscard]] bool isScanNumberValid(int scanNumber) const;


private:
    Err parse(QXmlStreamReader &reader);
    Err parseMsRun(QXmlStreamReader &reader);
    Err parseScan(QXmlStreamReader &reader);


public:
    QScopedPointer<QFile> m_file;
    QMap<ScanNumber, MsScanInfo> *m_msScanInfo;
    QMap<ScanNumber, ScanPoints> *m_scanPoints;

};

MsReaderMzML::PrivateData::PrivateData(
        QMap<ScanNumber, MsScanInfo> *msScanInfo,
        QMap<ScanNumber, ScanPoints> *scanPoints
        )
: m_msScanInfo(msScanInfo)
, m_scanPoints(scanPoints) {}

MsReaderMzML::PrivateData::~PrivateData() {
    closeFile();
}

Err MsReaderMzML::PrivateData::openFile(const QString &filename) {

    ERR_INIT

    e = closeFile(); ree;

    m_file.reset(new QFile(filename));
    e = ErrorUtils::isTrue(m_file->open(QIODevice::ReadOnly)); ree;

    QXmlStreamReader xmlReader(m_file.data());

    return parse(xmlReader);
}

Err MsReaderMzML::PrivateData::closeFile() {

    ERR_INIT

    m_msScanInfo->clear();
    m_scanPoints->clear();
    m_file.reset();

    ERR_RETURN
}

bool MsReaderMzML::PrivateData::isScanNumberValid(int scanNumber) const {
    return scanNumber >= 1 && scanNumber <= m_msScanInfo->size();
}

namespace {

    bool isCurrentTag(const QXmlStreamReader &reader, const QString &tagName) {
        return (0 == QStringRef::compare(reader.name(), tagName, Qt::CaseInsensitive));
    }

}//NAMESPACE
Err MsReaderMzML::PrivateData::parse(QXmlStreamReader &reader) {

    ERR_INIT

    while (!reader.atEnd() && !reader.hasError()) {

        reader.readNext();

        if (reader.isStartDocument()) {
            continue;
        }

        if (reader.isStartElement()) {
            if (isCurrentTag(reader, mzmlElementName)) {
                continue;
            } else if (isCurrentTag(reader, spectrumListElementName)) {
                e = parseMsRun(reader); ree;
            }
        }
    }

    reader.clear();

    ERR_RETURN
}

Err MsReaderMzML::PrivateData::parseMsRun(QXmlStreamReader &reader) {

    ERR_INIT

    const QString SCAN_COUNT_ATTRIBUTE_NAME = QStringLiteral("count");

    bool ok = true;

    QXmlStreamAttributes attributes = reader.attributes();
    const int scanCount = attributes.value(SCAN_COUNT_ATTRIBUTE_NAME).toInt(&ok);
    e = ErrorUtils::isTrue(ok); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Scan Count:" << scanCount;

    while (!reader.atEnd() && !reader.hasError()) {

        reader.readNext();

        if (reader.isStartElement()) {
            if (isCurrentTag(reader, spectrumElementName)) {
                e = parseScan(reader); ree;
                continue;
            }

        }

        if (reader.isEndElement()) {
            if (isCurrentTag(reader, mzmlElementName)) {
                break;
            }
        }
    }

    e = ErrorUtils::isEqual(scanCount, m_msScanInfo->size()); eee_absorb;

    ERR_RETURN
}

namespace {

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

}//namespace
Err MsReaderMzML::PrivateData::parseScan(QXmlStreamReader &reader) {

    ERR_INIT

    MsScanInfo msScanInfo;

    e = getScanNumber(reader, &msScanInfo.scanNumber); ree;

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
                e = ErrorUtils::isTrue(ok); ree;
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
                e = processScanData(reader, &msScanInfo);  ree;
            }

            if (stringMatch(currentReaderName, PRECURSOR)) {
                e = processPrecursor(reader, &msScanInfo); ree;
            }

            if (stringMatch(currentReaderName, BINARY_DATA_ARRAY)) {
                QVector<ScanPoint> scanPoints;
                e = processBinaryData(reader, &scanPoints); ree;
                (*m_scanPoints)[msScanInfo.scanNumber] = scanPoints;
            }

            continue;
        }

        if (reader.isEndElement()) {
            if (isCurrentTag(reader, spectrumElementName)) {

                m_msScanInfo->insert(msScanInfo.scanNumber, msScanInfo);

                if (msScanInfo.scanNumber % 1000 == 0) {
                    qDebug() << "mzML Scans Parsed" << msScanInfo.scanNumber;
                }

                ERR_RETURN
            }
        }
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

MsReaderMzML::MsReaderMzML() {
    m_d.reset(new PrivateData(
            &this->m_msScanInfo,
            &this->m_scanPoints)
            );
}

MsReaderMzML::~MsReaderMzML() {
}

Err MsReaderMzML::openFile(const QString &filePath) {

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

Err MsReaderMzML::closeFile() {
    return m_d->closeFile();
}

MsReaderBase MsReaderMzML::msReaderBase() {

    MsReaderBase msReaderBase;
    // NOTE: it is important to set scanInfos first or setScanPoints will Err out.
    msReaderBase.setMsScanInfo(*m_d->m_msScanInfo);
    msReaderBase.setScanPoints(*m_d->m_scanPoints);

    return msReaderBase;
}
