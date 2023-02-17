#include "PeptidesLibraryTron.h"
#include "Error.h"

#include <QtTest>

class PeptidesLibraryTronTests : public QObject
{
    Q_OBJECT

public:

    PeptidesLibraryTronTests() = default;
    ~PeptidesLibraryTronTests() override = default;

private slots:
    void processFastaTest();
    void removeDuplicatesFromPeptideTest();
    void addDecoysTest();
    void addVariableModificationsToPeptidesTest();

private:

    PythiaParameters pythiaParameters() {

        Modification modOx('M', "Ox", ModificationType::DYNAMIC, "O");
        Modification modDeam('N', "Deam", ModificationType::DYNAMIC, "N");
        Modification modCAM('C', "CAM", ModificationType::FIXED, "C2H3NO");

        PythiaParameters pythiaParameters;
        pythiaParameters.modifications.append({modCAM, modOx, modDeam});
        pythiaParameters.cTermCleavePoints.append({"K", "R"});
        pythiaParameters.allowedMissedCleavages = 1;
        pythiaParameters.addDecoys = true;
        pythiaParameters.maxModificationsPeptide = 2;

        return pythiaParameters;
    }

    QString fastaFilePath() {
        return QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");
    }
};

void PeptidesLibraryTronTests::processFastaTest()
{
    ERR_INIT;

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);

    QCOMPARE(pepLib.m_peptideSequences.size(), 87612);
}

void PeptidesLibraryTronTests::removeDuplicatesFromPeptideTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 87612);

    e = pepLib.removeDuplicatesFromPeptide();
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 85574);
}

void PeptidesLibraryTronTests::addDecoysTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 87612);

    e = pepLib.removeDuplicatesFromPeptide();
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 85574);

    e = pepLib.addDecoys(666);
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 171126);
}

void PeptidesLibraryTronTests::addVariableModificationsToPeptidesTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);
    QCOMPARE(pepLib.m_peptideSequences.size(), 87612);

    pepLib.addVariableModificationsToPeptides();


}


QTEST_MAIN(PeptidesLibraryTronTests)

#include "PeptidesLibraryTronTests.moc"
