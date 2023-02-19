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

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    FragLibraryReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

}

void FragLibraryReaderTests::writeFileTest() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    FragLibraryReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

    const QString writeFileName = QStringLiteral("SoLetItBeWritten.parquet");
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
