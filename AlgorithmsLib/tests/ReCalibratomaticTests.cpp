#include "ReCalibratomatic.h"
#include "GlobalSettings.h"
#include "Error.h"

#include <QtTest>

#include <iostream>

class ReCalibratomaticTests : public QObject
{
    Q_OBJECT

public:

    ReCalibratomaticTests() = default;
    ~ReCalibratomaticTests() override = default;

private slots:

    void svmTest();

private:

    static QString fragFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("test" + S_GLOBAL_SETTINGS.DOT_FRAGLIB);
    }
};


void ReCalibratomaticTests::svmTest() {

    ReCalibratomatic::svmTest();

}


QTEST_MAIN(ReCalibratomaticTests)

#include "ReCalibratomaticTests.moc"
