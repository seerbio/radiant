#include <QtTest>
#include <QCoreApplication>
#include <QDebug>

#include "Error.h"
#include "MsScansDenoiseTron.h"


using namespace Error;


class MsScansDenoiseTronTests : public QObject
{
    Q_OBJECT

public:

    MsScansDenoiseTronTests() = default;
    ~MsScansDenoiseTronTests() override = default;

private Q_SLOTS:

    void denoiseMsScans();

};

void MsScansDenoiseTronTests::denoiseMsScans() {



}


QTEST_MAIN(MsScansDenoiseTronTests)

#include "MsScansDenoiseTronTests.moc"
