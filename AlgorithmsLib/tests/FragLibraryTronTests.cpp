#include "FragLibraryTron.h"
#include "Error.h"
#include "PythiaParameterReader.h"

#include <QtTest>

#include <iostream>

class FragLibraryTronTests : public QObject
{
    Q_OBJECT

public:

    FragLibraryTronTests() = default;
    ~FragLibraryTronTests() override = default;

private slots:

    void readWriteTest();


private:

    QVector<FragLibIon> buildTestData() {

        std::srand(666);
        const int n = 10000;

        std::mt19937 re(std::random_device {}());

        std::uniform_real_distribution<double> unifMass(1000.0, 4800.0);
        auto generatorMass = bind(unifMass, std::ref(re));

        QVector<double> randoMass(n);
        std::generate(randoMass.begin(), randoMass.end(), std::ref(generatorMass));

        std::uniform_real_distribution<double> unifMz(400.0, 1800.0);
        auto generatorMz = bind(unifMz, std::ref(re));

        QVector<double> randoMz(n);
        std::generate(randoMz.begin(), randoMz.end(), std::ref(generatorMz));

        QVector<FragLibIon> fragLibIons;
        for (int i = 0; i < n; i++) {
            fragLibIons.push_back({1, randoMz.at(i), randoMass.at(i)});
        }

        return fragLibIons;
    }


private:

    static QString fragFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("test.fragLib");
    }
};

void FragLibraryTronTests::readWriteTest()
{

    ERR_INIT;

    const QVector<FragLibIon> fragLibIons = buildTestData();

    e = FragLibraryTron::writeFragLibIons(
            fragLibIons,
            fragFilePath()
            );
    QCOMPARE(e, eNoError);

    FragLibraryTron fragLibTron;
    e = fragLibTron.readFragLibIons(fragFilePath());
    QCOMPARE(e, eNoError);

    FragLibIon previous = FragLibIon();
    for (const FragLibIon &fl : fragLibTron.m_fragLibIons) {
        QVERIFY(previous.mzFrag <= fl.mzFrag);
        previous = fl;
    }

    QFile::remove(fragFilePath());
    QFileInfo checkFile2(fragFilePath());
    QCOMPARE(checkFile2.exists(), false);
}


QTEST_MAIN(FragLibraryTronTests)

#include "FragLibraryTronTests.moc"
