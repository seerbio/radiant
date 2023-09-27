#include "FDRCLassifierNeuralNet.h"

#include <QtTest/QtTest>

#include <iostream>

class FDRCLassifierNeuralNetTests : public QObject
{
    Q_OBJECT

public:
    FDRCLassifierNeuralNetTests() = default;
    ~FDRCLassifierNeuralNetTests() override = default;

private Q_SLOTS:
    void initTest();
};

void FDRCLassifierNeuralNetTests::initTest() {

    ERR_INIT

    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = fdrClassifierNeuralNet.init(
            PythiaParameters(),
            1,
            5,
            0.1,
            0.01
            );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            0,
            5,
            0.1,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            1,
            0,
            0.1,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            1,
            1,
            0.0,
            0.01
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            1,
            1,
            0.1,
            0.0
    );
    QCOMPARE(e, eError);

    e = fdrClassifierNeuralNet.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            1,
            1,
            0.1,
            0.0000001
    );
    QCOMPARE(e, eNoError);

}






QTEST_MAIN(FDRCLassifierNeuralNetTests)
#include "FDRCLassifierNeuralNetTests.moc"
