#include "BiophysicalCalcs.h"

#include "AminoAcids.h"

#include <QHash>
#include <QtTest>

class BiophysicalCalcsTests : public QObject
{
    Q_OBJECT

public:
    BiophysicalCalcsTests() = default;
    ~BiophysicalCalcsTests() override = default;

private slots:

    void calculatePeptideMassTest();
    void calculateMassFromThomsonTest();


};


void BiophysicalCalcsTests::calculatePeptideMassTest() {

    const QString seq = QStringLiteral("TT");
    const QHash<int, double> mod = {
            {0, 100.0}
    };

    const double seqMass = BiophysicalCalcs::calculatePeptideMass(
            seq,
            AminoAcids(),
            mod
            );

    const double expectedSeqMass = 320.106;

    const auto test = static_cast<int>(seqMass);
    const auto expected = static_cast<int>(expectedSeqMass);

    QCOMPARE(test, expected);
}

void BiophysicalCalcsTests::calculateMassFromThomsonTest() {

    const double result = BiophysicalCalcs::calculateMassFromThomson(501.5, 2, 1);
    QCOMPARE(result, 999.98234706);

}


QTEST_MAIN(BiophysicalCalcsTests)

#include "BiophysicalCalcsTests.moc"
