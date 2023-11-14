#include "KarnnNeuralNet.h"

#include "MathUtils.h"

#include <QElapsedTimer>
#include <QtTest/QtTest>

#include <iostream>
#include <fstream>


bool readDataCSV(
        const std::string &csvFilepath,
        QVector<QVector<double>> *dataVec,
        QStringList *seqsList,
        QVector<bool> *isDecoy,
        int *rows,
        int *cols
) {

    std::ifstream inputFile(csvFilepath);

    if (!inputFile.is_open()) {
        qDebug() << "Unable to open file for reading.";
        return false;
    }

    std::string line;
    int rowCounter = 0;
    while (std::getline(inputFile, line)) {

        QStringList stringSplit = QString::fromStdString(line).split(",", Qt::SkipEmptyParts);

        const QString seq = stringSplit.front();
        stringSplit.pop_front();
        seqsList->push_back(seq);

        bool ok;
        const float result = stringSplit.back().toFloat(&ok);
        stringSplit.pop_back();
        if (!ok) {
            throw;
        }
        if (result == 1) {
            isDecoy->push_back(false);
        }
        else {
            isDecoy->push_back(true);
        }

        *cols = stringSplit.size();

        QVector<double> v;
        for (const QString &val : stringSplit) {
            v.push_back(val.toDouble(&ok));
            if (!ok) {
                throw;
            }
        }

        dataVec->push_back(v);

        rowCounter++;
    }

    inputFile.close();
    *rows = rowCounter;

    return true;

}


class KarnnNeuralNetTests : public QObject
{
    Q_OBJECT

public:
    KarnnNeuralNetTests() = default;
    ~KarnnNeuralNetTests() = default;


private Q_SLOTS:
    void initTest();
    void runTest();

};


void KarnnNeuralNetTests::initTest() {

    ERR_INIT

    KarnnNeuralNet karnnNeuralNet;

    e = karnnNeuralNet.init(
            0,
            5,
            0.003
            );
    QCOMPARE(e, eError);

    e = karnnNeuralNet.init(
            12,
            0,
            0.003
    );
    QCOMPARE(e, eError);

    e = karnnNeuralNet.init(
            12,
            5,
            0
    );
    QCOMPARE(e, eError);

    e = karnnNeuralNet.init(
            12,
            5,
            0.003
    );
    QCOMPARE(e, eNoError);

}

void KarnnNeuralNetTests::runTest() {

    ERR_INIT

    //TODO use proper paths and file storage
    const std::string dataPath = "/home/anichols/Desktop/Data/Reports/vec1.csv";

    QVector<QVector<double>> dataVec;
    QStringList seqsList;
    int rows;
    int cols;
    QVector<bool> isDecoy;
    const bool ranOK = readDataCSV(
            dataPath,
            &dataVec,
            &seqsList,
            &isDecoy,
            &rows,
            &cols
    );

    KarnnNeuralNet karnnNeuralNet;
    e = karnnNeuralNet.init(
            12,
            5,
            0.003
    );
    QCOMPARE(e, eNoError);

    const int epochs = 1;
    std::vector<float> nnScores;
    e = karnnNeuralNet.run(
            dataVec,
            isDecoy,
            epochs,
            &nnScores
            );
    QCOMPARE(e, eNoError);



}


QTEST_MAIN(KarnnNeuralNetTests)
#include "KarnnNeuralNetTests.moc"



