//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FragLibReaderRow.h"
#include "SpecLibReader.h"
#include "SpecLibSrc/SpecLibStructures.h"

#include <QtTest/QtTest>

#include "FragLibReader.h"

#include <fstream>

namespace {

    void writeString(std::ofstream &out, const std::string &value) {
        const int size = static_cast<int>(value.size());
        out.write(reinterpret_cast<const char*>(&size), sizeof(int));
        if (size > 0) {
            out.write(value.data(), size);
        }
    }

    template <typename T>
    void writeVector(std::ofstream &out, const std::vector<T> &values) {
        const int size = static_cast<int>(values.size());
        out.write(reinterpret_cast<const char*>(&size), sizeof(int));
        if (size > 0) {
            out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(size * sizeof(T)));
        }
    }

    void writeMinimalEntry(
        std::ofstream &out,
        int index,
        int charge,
        int length,
        float precursorMz
        ) {
        const float iRt = 1.0f;
        const float sRt = 0.1f;
        const float libQValue = 0.0f;
        const float iIM = 0.0f;
        const float sIM = 0.0f;
        out.write(reinterpret_cast<const char*>(&index), sizeof(int));
        out.write(reinterpret_cast<const char*>(&charge), sizeof(int));
        out.write(reinterpret_cast<const char*>(&length), sizeof(int));
        out.write(reinterpret_cast<const char*>(&precursorMz), sizeof(float));
        out.write(reinterpret_cast<const char*>(&iRt), sizeof(float));
        out.write(reinterpret_cast<const char*>(&sRt), sizeof(float));
        out.write(reinterpret_cast<const char*>(&libQValue), sizeof(float));
        out.write(reinterpret_cast<const char*>(&iIM), sizeof(float));
        out.write(reinterpret_cast<const char*>(&sIM), sizeof(float));

        Product fragment;
        fragment.mz = precursorMz / static_cast<float>(charge);
        fragment.height = 1.0f;
        fragment.charge = 1;
        fragment.type = 1;
        fragment.index = 1;
        fragment.loss = 0;
        writeVector(out, std::vector<Product>{fragment});

        const int hasDecoy = 0;
        const int entryFlags = 0;
        const int proteotypic = 0;
        const int pidIndex = 0;
        const float pgQValue = 0.0f;
        const float ptmQValue = 1.0f;
        const float siteConf = 1.0f;
        out.write(reinterpret_cast<const char*>(&hasDecoy), sizeof(int));
        out.write(reinterpret_cast<const char*>(&entryFlags), sizeof(int));
        out.write(reinterpret_cast<const char*>(&proteotypic), sizeof(int));
        out.write(reinterpret_cast<const char*>(&pidIndex), sizeof(int));
        writeString(out, "");
        out.write(reinterpret_cast<const char*>(&pgQValue), sizeof(float));
        out.write(reinterpret_cast<const char*>(&ptmQValue), sizeof(float));
        out.write(reinterpret_cast<const char*>(&siteConf), sizeof(float));
    }

    void writeSyntheticSpecLib(const QString &filePath) {
        std::ofstream out(filePath.toStdString(), std::ios::binary | std::ios::trunc);
        Q_ASSERT(out.good());

        const int version = -3;
        const int gd = 1;
        const int gc = 1;
        const int ip = 0;
        out.write(reinterpret_cast<const char*>(&version), sizeof(int));
        out.write(reinterpret_cast<const char*>(&gd), sizeof(int));
        out.write(reinterpret_cast<const char*>(&gc), sizeof(int));
        out.write(reinterpret_cast<const char*>(&ip), sizeof(int));

        writeString(out, "test");
        writeString(out, "test.fasta");

        const int proteinsSize = 0;
        out.write(reinterpret_cast<const char*>(&proteinsSize), sizeof(int));

        const int proteinGroupsSize = 1;
        out.write(reinterpret_cast<const char*>(&proteinGroupsSize), sizeof(int));
        const int sizeP = 0;
        out.write(reinterpret_cast<const char*>(&sizeP), sizeof(int));
        writeString(out, "");
        writeString(out, "PG1");
        writeString(out, "");
        writeVector(out, std::vector<int>{});
        writeVector(out, std::vector<int>{});
        writeVector(out, std::vector<int>{});

        const std::vector<std::string> precursors = {"ACDE2", "ACDX2"};
        const int precursorsSize = static_cast<int>(precursors.size());
        out.write(reinterpret_cast<const char*>(&precursorsSize), sizeof(int));
        for (const std::string &precursor : precursors) {
            writeString(out, precursor);
        }

        const int namesSize = 0;
        out.write(reinterpret_cast<const char*>(&namesSize), sizeof(int));

        const int genesSize = 0;
        out.write(reinterpret_cast<const char*>(&genesSize), sizeof(int));

        const double iRtMin = 0.0;
        const double iRtMax = 100.0;
        out.write(reinterpret_cast<const char*>(&iRtMin), sizeof(double));
        out.write(reinterpret_cast<const char*>(&iRtMax), sizeof(double));

        const int entriesSize = 2;
        out.write(reinterpret_cast<const char*>(&entriesSize), sizeof(int));
        writeMinimalEntry(out, 0, 2, 4, 500.0f);
        writeMinimalEntry(out, 1, 2, 4, 600.0f);

        const int elutionGroupsSize = 0;
        out.write(reinterpret_cast<const char*>(&elutionGroupsSize), sizeof(int));
    }

}

class SpecLibReaderTests : public QObject
{
    Q_OBJECT

public:
    SpecLibReaderTests() = default;
    ~SpecLibReaderTests() override = default;

private Q_SLOTS:

    static void getFragLibReaerRowsTest();
    static void peptideSearchTroubleShoot();
    static void getFragLibReaerRowsFiltersInvalidSequencesTest();

};


void SpecLibReaderTests::peptideSearchTroubleShoot() {

    QSKIP("Turned on only during troubleshooting");

    ERR_INIT

    const QString filePath = "/home/andrewnichols/Downloads/uniprot_homosapiens_isoforms_con_20221209.predicted(1).speclib";

    QList<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        filePath,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    const QString peptideSequence = "EAQEFVK";
    const int charge = 1;

    for (const FragLibReaderRow &flrr : fragLibReaderRows) {

        if (flrr.peptideSequenceChargeKey != peptideSequence + "|" + QString::number(charge)) {
            continue;
        }

        const QStringList ionLabels = flrr.ionLabels.split(S_GLOBAL_SETTINGS.SEPARATOR);
        for (int i = 0; i < flrr.mzVals.size(); i++) {
            qDebug() << flrr.mzVals.at(i) << flrr.intensityVals.at(i) << ionLabels.at(i);
        }
        qDebug() << "Precursor charge" << charge;
    }

}

void SpecLibReaderTests::getFragLibReaerRowsTest() {

    QSKIP("Build a test for this when a small lib file is found");

    ERR_INIT

    // const QString specLibFile = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/lib.predicted.speclib");
    const QString specLibFile
        = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.speclib");
    QList<FragLibReaderRow> fragLibReaderRows;
    e = SpecLibReader::getFragLibReaerRows(
        specLibFile,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    const QString fragLibFFFile
    = QStringLiteral("/home/andrewnichols/Desktop/Data/Libraries/diannformat-human_plasma_arath_entrapment-lib.tsv.mods.fragLibFF");
    QList<FragLibReaderRow> fragLibReaderRowsFF;
    e = FragLibReader::getFragLibReaderRows(
        fragLibFFFile,
        false,
        &fragLibReaderRowsFF
        );
    QCOMPARE(e, eNoError);

    qDebug() << fragLibReaderRowsFF.size() << fragLibReaderRows.size() << "Sizes";

    for (int i = 0; i < fragLibReaderRows.size(); i++) {
        QCOMPARE(fragLibReaderRows.at(i).peptideSequenceChargeKey, fragLibReaderRowsFF.at(i).peptideSequenceChargeKey);
        QCOMPARE(fragLibReaderRows.at(i).intensityVals, fragLibReaderRowsFF.at(i).intensityVals);
        // if (fragLibReaderRows.at(i).ionLabels != fragLibReaderRowsFF.at(i).ionLabels) {
        //     qDebug() << fragLibReaderRows.at(i).ionLabels << fragLibReaderRowsFF.at(i).ionLabels;
        //     qDebug() << fragLibReaderRows.at(i).intensityVals << fragLibReaderRowsFF.at(i).intensityVals;
        //
        // }


    }


}

void SpecLibReaderTests::getFragLibReaerRowsFiltersInvalidSequencesTest() {
    ERR_INIT

    const QString specLibFile = QDir(qApp->applicationDirPath()).filePath("SpecLibReaderTests.invalid_sequences.speclib");
    writeSyntheticSpecLib(specLibFile);

    QList<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
        specLibFile,
        false,
        &fragLibReaderRows
        );
    QCOMPARE(e, eNoError);

    QCOMPARE(fragLibReaderRows.size(), 1);
    QCOMPARE(fragLibReaderRows.front().peptideSequenceChargeKey, QString("ACDE|2"));
}


QTEST_MAIN(SpecLibReaderTests)
#include "SpecLibReaderTests.moc"
