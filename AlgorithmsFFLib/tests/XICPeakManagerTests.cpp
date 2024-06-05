//
// Created by anichols on 11/07/2021.
//

#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FragLibReader.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"

#include <QtTest/QtTest>

#include "MsCalibratomatic.h"

class XICPeakManagerTests : public QObject
{
    Q_OBJECT

public:
    XICPeakManagerTests() = default;
    ~XICPeakManagerTests() override = default;

private Q_SLOTS:

    static void findPeaksTest();


};

void XICPeakManagerTests::findPeaksTest() {

    QSKIP("This is for dev purposes.  Write real tests for this class");

    //TODO Write tests
}


QTEST_MAIN(XICPeakManagerTests)
#include "XICPeakManagerTests.moc"
