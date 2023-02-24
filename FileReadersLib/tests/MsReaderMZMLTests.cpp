//
// Created by anichols on 11/07/2021.
//

#include "MsReaderMzML.h"
#include "FastaReader.h"
#include "GlobalSettings.h"

#include <QString>
#include <QtTest/QtTest>
#include <QXmlStreamReader>



class MsReaderMZMLTests : public QObject
{
    Q_OBJECT

public:
    MsReaderMZMLTests() = default;
    ~MsReaderMZMLTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    void openFileTest();

private:

    //TODO use proper path procedures.
    const QString m_filepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

};


void MsReaderMZMLTests::openFileTest() {

//    QSKIP("Waiting for small file");
    ERR_INIT

    const QString cacheName = m_filepath + S_GLOBAL_SETTINGS.DOT_CACHE;

    MsReaderMzML reader;
    e = reader.openFile(m_filepath);
    QCOMPARE(e, Error::eNoError);

////    e = reader.createTandemScanIonsCache(cacheName);
////    QCOMPARE(e, Error::eNoError);
//
//    e = reader.readFromCache(cacheName);
//    QCOMPARE(e, Error::eNoError);
//
////    e = reader.buildUniqueTandemScanIons();
////    QCOMPARE(e, Error::eNoError);
//
//    qDebug() << "IONS SIZE" << reader.m_tandemScanIons.size();
//    qDebug() << "UNIQUES IONS SIZE" << reader.m_uniqueTandemScanIons.size();

}



QTEST_MAIN(MsReaderMZMLTests)
#include "MsReaderMZMLTests.moc"
