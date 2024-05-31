//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "MsReaderBase.h"

#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

#include <QtTest/QtTest>
#include <QtConcurrent/QtConcurrent>

#include "GlobalSettings.h"
#include "MsReaderBase.h"


class PlayGroundTests : public QObject
{
    Q_OBJECT

public:
    PlayGroundTests() = default;
    ~PlayGroundTests() override = default;

private Q_SLOTS:

void testme();

};

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

namespace {

    struct FileChunk {
        const char* data;
        size_t start;
        size_t end;
        std::vector<std::string> *strs;
    };

    bool isCurrentTag(const QXmlStreamReader &reader, const QString &tagName) {
        return (0 == QStringRef::compare(reader.name(), tagName, Qt::CaseInsensitive));
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

    std::string parseScan(QXmlStreamReader &reader, const FileChunk &fc) {

        ERR_INIT

        MsScanInfo msScanInfo;

        e = getScanNumber(reader, &msScanInfo.scanNumber);
        return QString::number(msScanInfo.scanNumber).toStdString();

    //     while (!reader.atEnd() && !reader.hasError()) {
    //
    //         reader.readNext();
    //
    //         if (reader.isStartElement()) {
    //
    //             const QXmlStreamAttributes attributesCurrentTag = reader.attributes();
    //             const QStringRef currentAccession = attributesCurrentTag.value(ACCESSION);
    //             const QStringRef currentValue = attributesCurrentTag.value(VALUE);
    //             const QStringRef currentReaderName = reader.name();
    //
    //             bool ok = true;
    //
    //             if (stringMatch(currentAccession, MS_LEVEL)) {
    //                 msScanInfo.msLevel = currentValue.toInt(&ok);
    //                 e = ErrorUtils::isTrue(ok); ree;
    //             }
    //
    // //            if (stringMatch(currentAccession, FAIMS_VOLTAGE)) {
    // //                msScanInfo.faimsVoltage = static_cast<int>(currentValue.toDouble(&ok));
    // //                e = ErrorUtils::isTrue(ok); ree;
    // //            }
    //
    // //            if (stringMatch(currentAccession, TIC)) {
    // //                msScanInfo.TIC = currentValue.toDouble(&ok);
    // //                e = ErrorUtils::isTrue(ok); ree;
    // //            }
    // //
    // //            if (stringMatch(currentAccession, BASEPEAK_INTENSITY)) {
    // //                msScanInfo.basePeakIntensity = currentValue.toDouble(&ok);
    // //                e = ErrorUtils::isTrue(ok); ree;
    // //            }
    //
    //             // if (stringMatch(currentReaderName, SCAN)) {
    //             //     e = processScanData(reader, &msScanInfo);  ree;
    //             // }
    //             //
    //             // if (stringMatch(currentReaderName, PRECURSOR)) {
    //             //     e = processPrecursor(reader, &msScanInfo); ree;
    //             // }
    //             //
    //             // if (stringMatch(currentReaderName, BINARY_DATA_ARRAY)) {
    //             //     QVector<ScanPoint> scanPoints;
    //             //     e = processBinaryData(reader, &scanPoints); ree;
    //             //     (*m_scanPoints)[msScanInfo.scanNumber] = scanPoints;
    //             // }
    //
    //             continue;
    //         }
    //
    //         if (reader.isEndElement()) {
    //             if (isCurrentTag(reader, spectrumElementName)) {
    //
    //                 // m_msScanInfo->insert(msScanInfo.scanNumber, msScanInfo);
    //
    //                 if (msScanInfo.scanNumber % 1000 == 0) {
    //                     qDebug() << "mzML Scans Parsed" << msScanInfo.scanNumber;
    //                 }
    //
    //                 ERR_RETURN
    //             }
    //         }
    //     }
    //
    //     ERR_RETURN
    }

    std::vector<std::string> parseMsRun(QXmlStreamReader &reader, const FileChunk& chunk) {

        ERR_INIT

        const QString SCAN_COUNT_ATTRIBUTE_NAME = QStringLiteral("count");

        bool ok = true;

        QXmlStreamAttributes attributes = reader.attributes();
        const int scanCount = attributes.value(SCAN_COUNT_ATTRIBUTE_NAME).toInt(&ok);
        e = ErrorUtils::isTrue(ok);

        qDebug() << "Scan Count:" << scanCount;
        std::vector<std::string> scanNumbers;
        while (!reader.atEnd() && !reader.hasError()) {

            reader.readNext();

            if (reader.isStartElement()) {
                if (isCurrentTag(reader, spectrumElementName)) {
                     scanNumbers.push_back(parseScan(reader, chunk));
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

        return scanNumbers;
    }

    std::vector<std::string> processChunk(const FileChunk& chunk) {
        size_t elementCount = 0;
        const size_t maxSegmentSize = 1024 * 1024;  // 1 MB segments

        std::vector<std::string> allstring;

        try {
            for (size_t i = chunk.start; i < chunk.end; i += maxSegmentSize) {
                size_t segmentEnd = std::min(i + maxSegmentSize, chunk.end);
                QString segmentContent = QString::fromUtf8(chunk.data + i, segmentEnd - i);
                QXmlStreamReader xmlReader(segmentContent);

                while (!xmlReader.atEnd()) {

                    xmlReader.readNext();

                    if (xmlReader.isStartDocument()) {
                        continue;
                    }

                    if (xmlReader.isStartElement()) {
                        if (isCurrentTag(xmlReader, mzmlElementName)) {
                            continue;
                        } else if (isCurrentTag(xmlReader, spectrumListElementName)) {
                            std::vector<std::string> asdf = parseMsRun(xmlReader, chunk); //ree;
                            for (auto &ss : asdf) {
                                allstring.push_back(ss);
                            }
                        }
                    }

                }

                if (xmlReader.hasError() && !xmlReader.atEnd()) {
                    qWarning() << "XML error:" << xmlReader.errorString();
                }
            }
        } catch (const std::exception& e) {
            qCritical() << "Exception caught in processChunk:" << e.what();
            throw;
        }

        return allstring;

    }

    std::vector<std::string> processChunkCount(const FileChunk& chunk) {

        std::vector<std::string> ids;

        std::string str;
        for (size_t i = chunk.start; i < chunk.end; ++i) {

            str += chunk.data[i];

            if (chunk.data[i] == '\n') {

                if (QString::fromStdString(str).contains("<binary>")) {
                    ids.push_back(str);
                }

                str.clear();
            }
        }

        return ids;
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

    void processFileInChunks(const QString& filename) {
        QFile file(filename);
        if (!file.open(QIODevice::ReadOnly)) {
            qCritical() << "Error opening file:" << file.errorString();
            return;
        }

        QFileInfo fileInfo(file);
        size_t fileSize = fileInfo.size();
        uchar* ucharData = file.map(0, fileSize);
        if (!ucharData) {
            qCritical() << "Error mapping file:" << file.errorString();
            return;
        }
        const char* data = reinterpret_cast<const char*>(ucharData);

        std::vector<std::string> allstring;

        const int numChunks = 8;
        size_t chunkSize = fileSize / numChunks;
        QVector<FileChunk> chunks;

        for (int i = 0; i < numChunks; ++i) {
            size_t start = i * chunkSize;
            size_t end = (i == numChunks - 1) ? fileSize : (i + 1) * chunkSize;
            // Adjust end to avoid splitting a line
            if (i != numChunks - 1) {
                while (end < fileSize && data[end] != '\n') {
                    ++end;
                }
            }
            chunks.append({data, start, end, &allstring});
        }

        // for (const FileChunk &fc : chunks) {
        //     std::vector<std::string> strs = processChunkCount(fc);
        //     for (const auto &s : strs) {
        //
        //         std::cout << s << std::endl;
        //     }
        //
        // }

        QFuture<std::vector<std::string>> future = QtConcurrent::mapped(chunks, processChunkCount);
        future.waitForFinished();

        qDebug() << allstring.size();

        int count = 0;
        for (const std::vector<std::string> &ss : future) {
            for (const auto &s : ss) {
                QString qs = QString::fromStdString(s);
                qs = qs.replace("\n", "");

                // QMap<QString, QString> res = parseAttributes(qs);
                // std::cout << qs.toStdString() << std::endl;
                // std::cout << res.value("index").toStdString() << std::endl;
                count++;
            }
        }
        std::cout << count << std::endl;

        file.unmap(ucharData);
    }

}

void PlayGroundTests::testme() {

    const char* filename  = "/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
    // const char* filename = "/home/anichols/Repos/PythiaDIACpp/FileReadersLib/tests/TestFiles/1min.mzML";
    int numThreads = 8;
    processFileInChunks(filename);

}

QTEST_MAIN(PlayGroundTests)
#include "PlayGroundTests.moc"
