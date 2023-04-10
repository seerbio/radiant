#include "FragLibraryTron.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "PythiaParameterReader.h"

#include <QtTest>

#include <iostream>

class FragLibraryTronDIATests : public QObject
{
    Q_OBJECT

public:

    FragLibraryTronDIATests() = default;
    ~FragLibraryTronDIATests() override = default;

private slots:

    void readTest();


private:


private:

    static QString fragFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("test" + S_GLOBAL_SETTINGS.DOT_FRAGLIB);
    }
};

void FragLibraryTronDIATests::readTest()
{

    ERR_INIT;


}


QTEST_MAIN(FragLibraryTronDIATests)

#include "FragLibraryTronDIATests.moc"
