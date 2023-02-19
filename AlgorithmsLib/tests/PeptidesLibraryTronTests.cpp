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
    void addTerminalModificationsToPeptidesTest();
    void addPeptideIdToPeptidesTest();
    void addMassToPeptidesTest();
    void execTest();
    void readWriteTest();

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

namespace {
    void reverseSortPeptidesBySeqWithMods(QVector<Peptide> *peptides) {
        std::sort(peptides->rbegin(),
                  peptides->rend(),
                  [](const Peptide &l, const Peptide &r){
                      return l.peptideStringWithMods() < r.peptideStringWithMods();
                  });
    }
}//namespace
void PeptidesLibraryTronTests::addVariableModificationsToPeptidesTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());

    Peptide pep;
    pep.sequence = "QMASSMSSNSLKSNDMSK";
    pep.previousResidue = '*';

    pepLib.m_peptides.clear();
    pepLib.m_peptides.push_back(pep);

    e = pepLib.addVariableModificationsToPeptides();
    QCOMPARE(e, eNoError);

    const QStringList expectedResults = {
            "QM[15.9949]ASSM[15.9949]SSNSLKSNDMSK",
            "QM[15.9949]ASSMSSN[14.0031]SLKSNDMSK",
            "QM[15.9949]ASSMSSNSLKSN[14.0031]DMSK",
            "QM[15.9949]ASSMSSNSLKSNDM[15.9949]SK",
            "QM[15.9949]ASSMSSNSLKSNDMSK",
            "QMASSM[15.9949]SSN[14.0031]SLKSNDMSK",
            "QMASSM[15.9949]SSNSLKSN[14.0031]DMSK",
            "QMASSM[15.9949]SSNSLKSNDM[15.9949]SK",
            "QMASSM[15.9949]SSNSLKSNDMSK",
            "QMASSMSSN[14.0031]SLKSN[14.0031]DMSK",
            "QMASSMSSN[14.0031]SLKSNDM[15.9949]SK",
            "QMASSMSSN[14.0031]SLKSNDMSK",
            "QMASSMSSNSLKSN[14.0031]DM[15.9949]SK",
            "QMASSMSSNSLKSN[14.0031]DMSK",
            "QMASSMSSNSLKSNDM[15.9949]SK",
            "QMASSMSSNSLKSNDMSK"
    };

    reverseSortPeptidesBySeqWithMods(&pepLib.m_peptides);
    QStringList results;
    for (const Peptide &p : pepLib.m_peptides) {
        results.push_back(p.peptideStringWithMods());
    }

    QCOMPARE(results, expectedResults);
}

void PeptidesLibraryTronTests::addTerminalModificationsToPeptidesTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());

    Peptide pep;
    pep.sequence = "QMASSMSSNSLKSNDMSK";
    pep.previousResidue = '*';

    Peptide pep2;
    pep2.sequence = "NSLKSNDMSK";
    pep2.previousResidue = 'S';

    pepLib.m_peptides.clear();
    pepLib.m_peptides.append({pep, pep2});

    e = pepLib.addTerminalModificationsToPeptides();
    QCOMPARE(e, eNoError);

    const QStringList expectedResults = {
            "Q[45.034]MASSMSSNSLKSNDMSK",
            "QMASSMSSNSLKSNDMSK",
            "NSLKSNDMSK"
    };

    reverseSortPeptidesBySeqWithMods(&pepLib.m_peptides);

    QStringList results;
    for (const Peptide &p : pepLib.m_peptides) {
        results.push_back(p.peptideStringWithMods());
    }

    QCOMPARE(results, expectedResults);
}

void PeptidesLibraryTronTests::addPeptideIdToPeptidesTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);

    e = pepLib.buildPeptides();
    QCOMPARE(e, eNoError);

    e = pepLib.addPeptideIdToPeptides();
    QCOMPARE(e, eNoError);
}

void PeptidesLibraryTronTests::addMassToPeptidesTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.processFasta(fastaFilePath());
    QCOMPARE(e, eNoError);

    e = pepLib.buildPeptides();
    QCOMPARE(e, eNoError);

    e = pepLib.addMassToPeptides();
    QCOMPARE(e, eNoError);

}

void PeptidesLibraryTronTests::execTest() {

    ERR_INIT

    PeptidesLibraryTron pepLib(pythiaParameters());
    e = pepLib.exec(fastaFilePath(), 666);
    QCOMPARE(e, eNoError);
}

void PeptidesLibraryTronTests::readWriteTest() {

    ERR_INIT

    QHash<ResidueIndex, ModificationMass> mods = {
            {1, 100.0},
            {2, 200.0}
    };

    Peptide pep1;
    pep1.id = 1;
    pep1.sequence = "SFSDDFS";
    pep1.mass = 1000;
    pep1.isDecoy = true;
    pep1.previousResidue = 'K';
    pep1.postResidue = 'Y';
    pep1.modifications = mods;

    QVector<Peptide> peptides;
    for (int i = 0; i < 100; i++){

        Peptide pep2 = pep1;
        pep2.id = i;
        peptides.append(pep2);
    }

    const QString libFilePath = fastaFilePath() + S_GLOBAL_SETTINGS.DOT_PEPLIB;

    e = PeptidesLibraryTron::writePeptidesLib(
            peptides,
            libFilePath
            );
    QCOMPARE(e, eNoError);

    PeptidesLibraryTron libTron;
    e = libTron.readPeptidesLib(libFilePath);
    QCOMPARE(e, eNoError);

    const Peptide &pep = libTron.m_peptidesLookupByPeptideId.value(66);
    QCOMPARE(pep.id, 66);
    QCOMPARE(pep.sequence, "SFSDDFS");
    QCOMPARE(pep.modifications.value(1), 100.0);
    QCOMPARE(pep.modifications.value(2), 200.0);

    peptides.push_back(peptides.first());
    e = PeptidesLibraryTron::writePeptidesLib(
            peptides,
            libFilePath
    );
    QCOMPARE(e, eNoError);

    PeptidesLibraryTron libTron2;
    e = libTron2.readPeptidesLib(libFilePath);
    QCOMPARE(e, eError);
}


QTEST_MAIN(PeptidesLibraryTronTests)

#include "PeptidesLibraryTronTests.moc"
