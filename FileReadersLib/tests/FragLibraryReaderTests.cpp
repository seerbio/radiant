//
// Created by anichols on 11/07/2021.
//

#include "FragLibraryReader.h"

#include <QtTest/QtTest>

class FragLibraryReaderTests : public QObject
{
    Q_OBJECT

public:
    FragLibraryReaderTests() = default;
    ~FragLibraryReaderTests() override = default;

private Q_SLOTS:

    void readFileTest();
    void writeFileTest();

};


void FragLibraryReaderTests::readFileTest() {

    const QString &fragLibIonsFilePath = QDir(qApp->applicationDirPath()).filePath("FragLibIons.parquet");

    FragLibraryReader reader;

    bool e;
    e = reader.readFile(fragLibIonsFilePath.toStdString());
    QCOMPARE(e, true);

    FragLibIon fliExpected;
    fliExpected.peptideId = 1;
    fliExpected.mzFrag = 66.6;
    fliExpected.peptideMass = 666.6;

    const FragLibIon &result = reader.m_fragLibIons.back();

    QCOMPARE(result.peptideId , fliExpected.peptideId);
    QCOMPARE(result.mzFrag , fliExpected.mzFrag);
    QCOMPARE(result.peptideMass , fliExpected.peptideMass);

}

void FragLibraryReaderTests::writeFileTest() {

    FragLibIon fli1;
    fli1.peptideId = 0;
    fli1.mzFrag = 66.6;
    fli1.peptideMass = 666.6;

    FragLibIon fli2;
    fli2.peptideId = 1;
    fli2.mzFrag = 66.6;
    fli2.peptideMass = 666.6;

    FragLibraryReader reader;
    reader.m_fragLibIons = {fli1, fli2};

    bool e;
    const QString writeFileName = QStringLiteral("SoLetItBeWrittenFragLibions.parquet");
    e = reader.writeFile(writeFileName.toStdString());
    QCOMPARE(e, true);

    QFileInfo checkFile(writeFileName);
    QCOMPARE(checkFile.exists(), true);

    QFile::remove(writeFileName);
    QFileInfo checkFile2(writeFileName);
    QCOMPARE(checkFile2.exists(), false);
}


QTEST_MAIN(FragLibraryReaderTests)
#include "FragLibraryReaderTests.moc"
