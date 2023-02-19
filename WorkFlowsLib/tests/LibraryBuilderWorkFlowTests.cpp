#include "LibraryBuilderWorkFlow.h"

#include "Error.h"
#include "PeptidesLibraryTron.h"
#include "PythiaParameterReader.h"

#include <QtTest>

class LibraryBuilderWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    LibraryBuilderWorkFlowTests() = default;
    ~LibraryBuilderWorkFlowTests() override = default;

private slots:

    void execTest();
    void testPeptideFragmentationTest();


private:

    static PythiaParameters pythiaParameters() {

        Modification modOx('M', "Ox", ModificationType::DYNAMIC, "O");
        Modification modDeam('N', "Deam", ModificationType::DYNAMIC, "N");
        Modification modCAM('C', "CAM", ModificationType::FIXED, "C2H3NO");
        Modification modAce("N-term-protein", "Ace", ModificationType::DYNAMIC, "C2H5O");

        PythiaParameters pythiaParameters;
        pythiaParameters.modifications.append({modCAM, modOx, modDeam, modAce});
        pythiaParameters.cTermCleavePoints.append({"K", "R"});
        pythiaParameters.allowedMissedCleavages = 1;
        pythiaParameters.addDecoys = true;
        pythiaParameters.maxModificationsPeptide = 2;

        PythiaParameterReader::applyFixedModificationsToAminoAcids(
                pythiaParameters,
                &pythiaParameters.aminoAcids
                );

        return pythiaParameters;
    }

    static QString fastaFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");
    }
};

void LibraryBuilderWorkFlowTests::execTest() {

    ERR_INIT

    LibraryBuilderWorkFlow libraryBuilderWorkFlow;

    e = libraryBuilderWorkFlow.exec(
            pythiaParameters(),
            fastaFilePath(),
            true
            );
    QCOMPARE(e, eNoError);

    const QString fragLibFilePath = fastaFilePath() + S_GLOBAL_SETTINGS.DOT_FRAGLIB;
    QFileInfo fi(fragLibFilePath);
    QCOMPARE(fi.exists(), true);

}

void LibraryBuilderWorkFlowTests::testPeptideFragmentationTest() {

    const QString peptideSeq = QStringLiteral("VTAK");
    QHash<ResidueIndex, ModificationMass> mods = {
            {1, 10.0}
    };

    const QVector<double> frags
        = LibraryBuilderWorkFlow::testPeptideFragmentation(peptideSeq, mods);

    const QVector<double> expectedFrags = {
            50.5415, 74.06, 100.076, 106.065, 109.579,
            141.584, 147.113, 165.102, 211.123, 218.15, 282.16, 329.198
    };

    QCOMPARE(frags.size(), expectedFrags.size());

    for (int i = 0; i < frags.size(); i++) {
        const int expectedMz = static_cast<int>(expectedFrags.at(i));
        const int foundMz = static_cast<int>(frags.at(i));
        QCOMPARE(foundMz, expectedMz);
    }

}


QTEST_MAIN(LibraryBuilderWorkFlowTests)

#include "LibraryBuilderWorkFlowTests.moc"
