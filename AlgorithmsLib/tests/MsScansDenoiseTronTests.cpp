#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "MsScansDenoiseTron.h"
#include "MsReaderMzML.h"


using namespace Error;


class MsScansDenoiseTronTests : public QObject
{
    Q_OBJECT

public:

    MsScansDenoiseTronTests() = default;
    ~MsScansDenoiseTronTests() override = default;

private Q_SLOTS:

    void denoiseMsScansTest();

private:

    //TODO use proper path procedures.
    const QString m_filepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");


};

void MsScansDenoiseTronTests::denoiseMsScansTest() {

    ERR_INIT

    MsReaderMzML reader;
    e = reader.openFile(m_filepath);
    QCOMPARE(e, Error::eNoError);



}


QTEST_MAIN(MsScansDenoiseTronTests)

#include "MsScansDenoiseTronTests.moc"
