#include "FeatureFinderHillBuilder.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderParquet.h"

#include <QtTest>

#include <iostream>

using namespace Error;


class FeatureFinderHillTests : public QObject
{
    Q_OBJECT

public:
    FeatureFinderHillTests() = default;
    ~FeatureFinderHillTests() override = default;

private slots:

    //TODO write tests


};


QTEST_MAIN(FeatureFinderHillTests)

#include "FeatureFinderHillTests.moc"
