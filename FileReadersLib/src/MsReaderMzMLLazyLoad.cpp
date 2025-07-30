
#include "MsReaderMzMLLazyLoad.h"

#include "ErrorUtils.h"
#include "ParallelUtils.h"
#include "SqlUtils.h"

#include <QFile>
#include <QMap>
#include <QXmlStreamReader>

#include <zlib.h>

#include <algorithm>
#include <iostream>
#include <QDebug>
#include <QFileInfo>

struct FileChunk {
    const char* data = nullptr;
    size_t start = -1;
    size_t end = -1;
    size_t overlap = -1;
};

const QString scanKey = QStringLiteral("scan");
const QString indexKey = QStringLiteral("index");

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


class Q_DECL_HIDDEN MsReaderMzMLLazyLoad::PrivateData  {

public:

    explicit PrivateData(QMap<ScanNumber, MsScanInfo> *msScanInfo);

    ~PrivateData();

    [[nodiscard]] Err openFile(const QString& filename) const;
    [[nodiscard]] Err closeFile() const;

    [[nodiscard]] bool isScanNumberValid(int scanNumber) const;

public:
    QMap<ScanNumber, MsScanInfo> *m_msScanInfo;

};

MsReaderMzMLLazyLoad::PrivateData::PrivateData(QMap<ScanNumber, MsScanInfo> *msScanInfo)
: m_msScanInfo(msScanInfo){}

MsReaderMzMLLazyLoad::PrivateData::~PrivateData() {
    Err e = closeFile();
}

namespace {

    Err memoryMapFile(
            QFile &file,
            qint64 *fileSize,
            uchar **ucharData
            ) {

        ERR_INIT

        if (!file.open(QIODevice::ReadOnly)) {
            qCritical() << "Error opening file:" << file.errorString();
            rrr(eFileError);
        }

        const QFileInfo fileInfo(file);
        *fileSize = fileInfo.size();
        *ucharData = file.map(0, *fileSize);
        if (!*ucharData) {
            qCritical() << "Error mapping file:" << file.errorString();
            rrr(eFileError);
        }

        ERR_RETURN
    }

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

		QString extractTagContents(
		std::ifstream& file,
		const std::streampos& tagPosition,
		const std::string& tag
		) {

		file.seekg(tagPosition);

		std::string line, tagContent;
		bool insideTag = false;
		std::string closingTag = "</" + tag.substr(1);

		while (std::getline(file, line)) {
			size_t start = (insideTag ? 0 : line.find(tag));
			if (start != std::string::npos) {

				start += tag.length();
				insideTag = true;
			}

			if (insideTag) {
				size_t close = line.find(closingTag);
				if (close != std::string::npos) {
					tagContent += line.substr(start, close - start);
					break;
				}

				tagContent += line.substr(start) + "\n";

			}
		}

		return QString::fromStdString(tagContent);
	}

	QPair<Err, std::streampos> findIndexListOffset(std::ifstream *file) {

		ERR_INIT

		const std::string& tag = "<indexListOffset>";

		std::streampos fileSize = file->tellg();
		std::streampos chunkSize = 4096;
		std::streampos position = fileSize;
		std::string buffer;

		std::string chunk;

		while (position > 0) {
			std::streampos readSize = (position >= chunkSize) ? chunkSize : position;
			position -= readSize;
			file->seekg(position);

			chunk.resize(readSize);
			file->read(&chunk[0], readSize);

			buffer.insert(0, chunk);

			size_t tagPos = buffer.rfind(tag);
			if (tagPos != std::string::npos) {

				const QString indexListOffsetLocationString
						= extractTagContents(*file, position + static_cast<std::streamoff>(tagPos), tag);

				return {e , static_cast<std::streamoff>(indexListOffsetLocationString.toLong())};
			}

			if (buffer.length() > tag.length()) {
				buffer = buffer.substr(0, tag.length());
			}
		}

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Error: Tag <indexListOffset> not found.";
		return {eFileError, {}};
	}

	Err initMsScanInfosWithOffsets(
		const QString &filePath,
		QMap<ScanNumber, MsScanInfo> *msScanInfos
		) {
		ERR_INIT

		e = ErrorUtils::fileExists(filePath); ree;

		std::ifstream file(filePath.toStdString(), std::ios::binary | std::ios::ate);
		e = ErrorUtils::isTrue(file.is_open()); ree;

		const QPair<Err, std::streampos> errVsStreamPos = findIndexListOffset(&file);
		e = errVsStreamPos.first; ree;
		const std::streampos posOffsets = errVsStreamPos.second;

		file.seekg(posOffsets);
		std::string line;

		while (std::getline(file, line)) {

			const QString lineQString = QString::fromStdString(line).trimmed();

			if (lineQString.mid(0, 7) != "<offset") {
				continue;
			}

			static QRegularExpression re(R"(<offset idRef="[^"]*">(\d+)</offset>)");
			QRegularExpressionMatch match = re.match(lineQString);

			e = ErrorUtils::isTrue(match.hasMatch()); ree;
			const QString offsetStr = match.captured(1);
			long offset;
			e = ErrorUtils::toLong(offsetStr, &offset); ree;

			static QRegularExpression re2(lineQString.contains(scanKey) ? R"(scan=(\d+))" : R"(index=(\d+))");
			QRegularExpressionMatch match2 = re2.match(lineQString);
			if (!match2.hasMatch()) {
				continue;
			}

			const QString scanNumberStr = match2.captured(1);

			int scanNumber;
			e = ErrorUtils::toInt(scanNumberStr, &scanNumber); ree;

			MsScanInfo msScanInfo;
			msScanInfo.scanOffsetStart = offset;
			msScanInfo.scanNumber = scanNumber;

			if (!msScanInfos->isEmpty()) {
				MsScanInfo *msScanInfoLast = &(*msScanInfos)[msScanInfos->lastKey()];
				msScanInfoLast->scanOffsetEnd = offset;
			}

			msScanInfos->insert(msScanInfo.scanNumber, msScanInfo);
		}

		ERR_RETURN
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
                    e = processSpectrumKey(
                    	str,
                    	&msScanInfoLocal.scanNumber
                    	); rree;
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
                    	const auto mzMinMax = std::minmax_element(mzValues.begin(), mzValues.end());
                    	msScanInfoLocal.mzMin = static_cast<float>(*mzMinMax.first);
                    	msScanInfoLocal.mzMax = static_cast<float>(*mzMinMax.second);
                    	msScanInfoLocal.pointCount = mzValues.size();
                    }
                    else if (spectrumType == SPECTRUM_TYPE::INTENSITY_ARR) {
                        intensityValues = decodedBLOB;
                    	const auto intensityValuesMinMax = std::minmax_element(intensityValues.begin(), intensityValues.end());
                    	msScanInfoLocal.intensityMin = static_cast<float>(*intensityValuesMinMax.first);
                    	msScanInfoLocal.intensityMax = static_cast<float>(*intensityValuesMinMax.second);
                    }
                }
                else if (str.contains(spectrumElementEndName)) {

                    if (mzValues.size() != intensityValues.size() || msScanInfoLocal.scanNumber < 0) {
                        msScanInfoLocal = MsScanInfo();
                        str.clear();
                        continue;
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
Err MsReaderMzMLLazyLoad::PrivateData::openFile(const QString &filename) const {

    ERR_INIT

    QFile file(filename);
    qint64 fileSize;
    uchar *ucharData;
    e = memoryMapFile(file, &fileSize, &ucharData); ree;

	e = initMsScanInfosWithOffsets(filename, m_msScanInfo); ree;
	m_msScanInfo->last().scanOffsetEnd = fileSize;

    const QVector<FileChunk> chunks = buildFileChunks(fileSize, ucharData);

#define RUN_PARALLEL
#ifdef RUN_PARALLEL
    QFuture<QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>>> future = QtConcurrent::mapped(chunks, processChunk);
    future.waitForFinished();

    for (const QPair<Err, QVector<QPair<MsScanInfo, ScanPoints>>> &result : future) {

    	e = result.first; ree;

    	for (const QPair<MsScanInfo, ScanPoints> &pr : result.second) {

        	const MsScanInfo &msScanInfo = pr.first;
            if (msScanInfo.scanNumber < 0) {
                continue;
            }

    		e = ErrorUtils::contains(msScanInfo.scanNumber, *m_msScanInfo); ree;

    		MsScanInfo *msScanInfoPntr = &(*m_msScanInfo)[msScanInfo.scanNumber];
    		msScanInfoPntr->frameIndex = msScanInfo.frameIndex;
    		msScanInfoPntr->collisionEnergy = msScanInfo.collisionEnergy;
    		msScanInfoPntr->intensityMin = msScanInfo.intensityMin;
    		msScanInfoPntr->intensityMax = msScanInfo.intensityMax;
    		msScanInfoPntr->pointCount = msScanInfo.pointCount;
    		msScanInfoPntr->msLevel = msScanInfo.msLevel;
    		msScanInfoPntr->scanTime = msScanInfo.scanTime;
    		msScanInfoPntr->precursorTargetMz = msScanInfo.precursorTargetMz;
    		msScanInfoPntr->isoWindowLower = msScanInfo.isoWindowLower;
    		msScanInfoPntr->isoWindowUpper = msScanInfo.isoWindowUpper;
    		msScanInfoPntr->ionMobilityDriftTime = msScanInfo.ionMobilityDriftTime;
    		msScanInfoPntr->mzMin = msScanInfo.mzMin;
    		msScanInfoPntr->mzMax = msScanInfo.mzMax;
    		msScanInfoPntr->ionMobilityIndex = msScanInfo.ionMobilityIndex;
        }
    }
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Scan Count" << m_msScanInfo->size();
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

Err MsReaderMzMLLazyLoad::PrivateData::closeFile() const {
    ERR_INIT
    m_msScanInfo->clear();
    ERR_RETURN
}

bool MsReaderMzMLLazyLoad::PrivateData::isScanNumberValid(int scanNumber) const {
    return scanNumber >= 1 && scanNumber <= m_msScanInfo->size();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

namespace {

    QPair<Err, QPair<MzVals, IntensityVals>> processChunkScanPoints(const FileChunk& chunk) {

        ERR_INIT

        ScanPoints scanPointsLocal;

        TYPES type = NULL_TYPE;
        COMPRESSION compression = COMPRESSION::NULL_COMPRESSION;
        SPECTRUM_TYPE spectrumType = SPECTRUM_TYPE::NULL_SPECTRUM_TYPE;

        MzVals mzValues;
        IntensityVals intensityValues;

    	const size_t length = chunk.end - chunk.start;
    	const QString qstr = QString::fromUtf8(chunk.data + chunk.start, length);
    	QStringList qstrList = qstr.split('\n', Qt::SkipEmptyParts);

        for (QString &str : qstrList) {

            str = str.trimmed();

        	// qDebug() << str << "SDFKJLS";

            if (str.contains(BIT64_FLOAT)) {
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

                QVector<float> decodedBLOB;
                if (type == TYPES::FLOAT64) {
                	const QVector<double> decodedBLOBDouble = SqlUtils::decodeBLOB<double>(binaryDataOut);
                	decodedBLOB.reserve(decodedBLOBDouble.size());
                	std::copy(decodedBLOBDouble.cbegin(), decodedBLOBDouble.cend(), std::back_inserter(decodedBLOB));
                }
                else if (type == TYPES::FLOAT32) {
                	decodedBLOB = SqlUtils::decodeBLOB<float>(binaryDataOut);
                }

                if (spectrumType == SPECTRUM_TYPE::MZ_ARR) {
                    mzValues = decodedBLOB;
                }
                else if (spectrumType == SPECTRUM_TYPE::INTENSITY_ARR) {
                    intensityValues = decodedBLOB;
                }
            }
            else if (str.contains(spectrumElementEndName)) {

                if (mzValues.size() != intensityValues.size()) {
                    scanPointsLocal = ScanPoints();
                    str.clear();
                    continue;
                }

                e = ErrorUtils::isEqual(mzValues.size(), intensityValues.size()); rree;
                scanPointsLocal.resize(mzValues.size());
                scanPointsLocal.reserve(mzValues.size());

                for (int j = 0; j < mzValues.size(); j++) {
                    const ScanPoint point(
                            static_cast<float>(mzValues.at(j)),
                            static_cast<float>(intensityValues.at(j))
                            );
                    scanPointsLocal[j] = point;
                }
                filterZeroIntensityQPoints(&scanPointsLocal);

                mzValues.resize(scanPointsLocal.size());
                intensityValues.resize(scanPointsLocal.size());
                mzValues.clear();
                intensityValues.clear();
                for (int j = 0; j < scanPointsLocal.size(); j++) {
                	const ScanPoint sp =scanPointsLocal[j];
                	mzValues[j] = sp.x();
                	intensityValues[j] = sp.y();
                }

                break;
            }

            str.clear();
        }

        return {e, {mzValues, intensityValues}};

    	// QXmlStreamReader xml(qstr);
    	// while (!xml.atEnd() && !xml.hasError()) {
	    //
    	// 	QXmlStreamReader::TokenType token = xml.readNext();
		   //  if (token == QXmlStreamReader::StartElement){
	    //
		   //  	if (xml.attributes().hasAttribute(ACCESSION)) {
		   //  		const QString accession = xml.attributes().value("accession").toString();
	    //
		   //  		if (accession == BIT64_FLOAT) {
		   //  			type = TYPES::FLOAT64;
		   //  		}
		   //  		else if (accession == BIT32_FLOAT) {
		   //  			type = TYPES::FLOAT32;
		   //  		}
		   //  		else if (accession == ZLIB_COMPRESSION) {
		   //  			compression = COMPRESSION::ZLIB;
		   //  		}
		   //  		else if (accession == NO_COMPRESSION) {
		   //  			compression = COMPRESSION::NO_COMPRESS;
		   //  		}
		   //  		else if (accession == MZ_ARRAY) {
		   //  			spectrumType = SPECTRUM_TYPE::MZ_ARR;
		   //  		}
		   //  		else if (accession == INTENSITY_ARRAY) {
		   //  			spectrumType = SPECTRUM_TYPE::INTENSITY_ARR;
		   //  		}
		   //  	}
	    //
    	// 		else if (xml.name() == "binary") {
	    //
    	// 			const QString binaryValue = xml.readElementText();
    	// 		    const std::string elementText = binaryValue.toStdString();
	    //
    	// 		    QByteArray binaryData = QByteArray::fromBase64(QByteArray::fromStdString(elementText));
	    //
    	// 		    QByteArray binaryDataOut;
    	// 		    if (compression == COMPRESSION::ZLIB) {
    	// 		        decompress(binaryData, &binaryDataOut);
    	// 		    }
    	// 		    else if (compression == COMPRESSION::NO_COMPRESS){
    	// 		        binaryDataOut = binaryData;
    	// 		    }
	    //
    	// 		    QVector<float> decodedBLOB;
    	// 		    if (type == TYPES::FLOAT64) {
    	// 			    const QVector<double> decodedBLOBDouble = SqlUtils::decodeBLOB<double>(binaryDataOut);
    	// 			    decodedBLOB.reserve(decodedBLOBDouble.size());
    	// 			    std::copy(decodedBLOBDouble.cbegin(), decodedBLOBDouble.cend(), std::back_inserter(decodedBLOB));
    	// 		    }
    	// 		    else if (type == TYPES::FLOAT32) {
    	// 			    decodedBLOB = SqlUtils::decodeBLOB<float>(binaryDataOut);
    	// 		    }
	    //
    	// 		    if (spectrumType == SPECTRUM_TYPE::MZ_ARR) {
    	// 		        mzValues = decodedBLOB;
    	// 		    }
    	// 		    else if (spectrumType == SPECTRUM_TYPE::INTENSITY_ARR) {
    	// 		        intensityValues = decodedBLOB;
    	// 		    }
    	// 		}
    	// 	}
    	// 	if (token == QXmlStreamReader::EndElement) {
    	// 		if (xml.name() == "spectrum") {
    	// 			if (mzValues.size() != intensityValues.size()) {
    	// 				scanPointsLocal = ScanPoints();
    	// 				continue;
    	// 			}
    	// 			e = ErrorUtils::isEqual(mzValues.size(), intensityValues.size()); rree;
	    //
    	// 			scanPointsLocal.resize(mzValues.size());
    	// 			scanPointsLocal.reserve(mzValues.size());
	    //
    	// 			for (int j = 0; j < mzValues.size(); j++) {
    	// 				const ScanPoint point(
					// 			static_cast<float>(mzValues.at(j)),
					// 			static_cast<float>(intensityValues.at(j))
					// 			);
    	// 				scanPointsLocal[j] = point;
    	// 			}
    	// 			filterZeroIntensityQPoints(&scanPointsLocal);
	    //
    	// 			mzValues.resize(scanPointsLocal.size());
    	// 			intensityValues.resize(scanPointsLocal.size());
    	// 			mzValues.clear();
    	// 			intensityValues.clear();
    	// 			for (int j = 0; j < scanPointsLocal.size(); j++) {
    	// 				const ScanPoint sp =scanPointsLocal[j];
    	// 				mzValues[j] = sp.x();
    	// 				intensityValues[j] = sp.y();
    	// 			}
	    //
    	// 			break;
    	// 		}
    	// 	}
    	// }
    }


}//namespace
Err MsReaderMzMLLazyLoad::extractScanPoints(
        const QString &fileName,
        const QVector<MsScanInfo*> &msScanInfos,
		QVector<MsScan> *msScans
        ) {
    ERR_INIT

    e = ErrorUtils::fileExists(fileName); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    QFile file(fileName);
    qint64 fileSize;
    uchar *ucharData;
    e = memoryMapFile(file, &fileSize, &ucharData); ree;

    const auto* data = reinterpret_cast<const char*>(ucharData);

    for (MsScanInfo *msi : msScanInfos) {
        size_t start = msi->scanOffsetStart;
        size_t end = msi->scanOffsetEnd;

        FileChunk fileChunk;
        fileChunk.data = data;
        fileChunk.start = start;
        fileChunk.end = end;
        fileChunk.overlap = 0;

        QPair<Err, QPair<MzVals, IntensityVals>> scanPointsResult = processChunkScanPoints(fileChunk);
        e = scanPointsResult.first; ree;

    	MsScan msScan;
    	msScan.msScanInfoPntr = msi;
    	msScan.mzVals = scanPointsResult.second.first;
    	msScan.intensityVals = scanPointsResult.second.second;

        msScans->push_back(msScan);
    }

	file.close();

    ERR_RETURN
}

// Err MsReaderMzMLLazyLoad::extractScanPoints(
// 		const QString &fileName,
// 		QVector<MsScanInfo> *msScanInfos,
// 		QVector<MsScan> *msScans
// 		) {
// 	ERR_INIT
//
// 	e = ErrorUtils::fileExists(fileName); ree;
// 	e = ErrorUtils::isNotEmpty(*msScanInfos); ree;
//
// 	QFile file(fileName);
// 	qint64 fileSize;
// 	uchar *ucharData;
// 	e = memoryMapFile(file, &fileSize, &ucharData); ree;
//
// 	const auto* data = reinterpret_cast<const char*>(ucharData);
//
// 	for (MsScanInfo &msi : *msScanInfos) {
// 		size_t start = msi.scanOffsetStart;
// 		size_t end = msi.scanOffsetEnd;
//
// 		FileChunk fileChunk;
// 		fileChunk.data = data;
// 		fileChunk.start = start;
// 		fileChunk.end = end;
// 		fileChunk.overlap = 0;
//
// 		QPair<Err, QPair<MzVals, IntensityVals>> scanPointsResult = processChunkScanPoints(fileChunk);
// 		e = scanPointsResult.first; ree;
//
// 		MsScan msScan;
// 		msScan.msScanInfoPntr = &msi;
// 		msScan.mzVals = scanPointsResult.second.first;
// 		msScan.intensityVals = scanPointsResult.second.second;
//
// 		msScans->push_back(msScan);
// 	}
//
// 	file.close();
//
// 	ERR_RETURN
// }

MsReaderMzMLLazyLoad::MsReaderMzMLLazyLoad() {
    m_d.reset(new PrivateData(&this->m_msScanInfo));
}

MsReaderMzMLLazyLoad::~MsReaderMzMLLazyLoad() {}

Err MsReaderMzMLLazyLoad::openFile(const QString &filePath) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(filePath); ree;
    m_filePath = filePath;

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            MzMLLazyLoadNamespace::MZML_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

    e = m_d->openFile(filePath); ree;

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    for (MsScanInfo &msi : m_msScanInfo) {
        const MzTargetKey mzTargetKey = msi.targetKey();
    	msi.frameIndex = m_mzTargetVsScanInfosPntrs[mzTargetKey].size();
        m_mzTargetVsScanInfosPntrs[mzTargetKey].push_back(&msi);

    	if (msi.msLevel == 1) {
    		m_mzMs1Min = std::min(m_mzMs1Min, msi.mzMin);
    		m_mzMs1Max = std::max(m_mzMs1Max, msi.mzMax);
    		m_intensityMs1Min = std::min(m_intensityMs1Min, msi.intensityMin);
    		m_intensityMs1Max = std::max(m_intensityMs1Max, msi.intensityMax);
    		continue;
    	}

    	m_mzMs2Min = std::min(m_mzMs2Min, msi.mzMin);
    	m_mzMs2Max = std::max(m_mzMs2Max, msi.mzMax);
    	m_intensityMs2Min = std::min(m_intensityMs2Min, msi.intensityMin);
    	m_intensityMs2Max = std::max(m_intensityMs2Max, msi.intensityMax);
    }

    e = printFileInfo(); ree;
    ERR_RETURN
}

Err MsReaderMzMLLazyLoad::closeFile() {
    return m_d->closeFile();
}

Err MsReaderMzMLLazyLoad::getMzTargetScanPoints(
    const MzTargetKey& targetKey,
    QVector<MsScan> *msScans
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetKey); ree;
    e = ErrorUtils::isTrue(isInit()); ree;

    QVector<MsScanInfo*> targetMsScanInfos;
    e = getMsScanInfos(targetKey, &targetMsScanInfos); ree;

    e = ErrorUtils::isNotEmpty(targetMsScanInfos); ree;

    e = extractScanPoints(
            targetMsScanInfos,
            msScans
            ); ree;

    ERR_RETURN
}

// MsReaderBase MsReaderMzMLLazyLoad::msReaderBase() const {
//
//     MsReaderBase msReaderBase;
//     msReaderBase.setMsScanInfo(*m_d->m_msScanInfo);
//     return msReaderBase;
// }

Err MsReaderMzMLLazyLoad::extractScanPoints(
	const QVector<MsScanInfo*> &msScanInfos,
	QVector<MsScan> *msScans
    ) {
    ERR_INIT
    e = MsReaderMzMLLazyLoad::extractScanPoints(
        m_filePath,
        msScanInfos,
        msScans
        ); ree;
    ERR_RETURN
}
