#include "ProteinDigestomatic.h"
#include "Error.h"
#include "FastaReader.h"

#include <QString>
#include <QStringList>
#include <QtTest>

class ProteinDigestomaticTests : public QObject
{
    Q_OBJECT

public:
    ProteinDigestomaticTests() = default;
    ~ProteinDigestomaticTests() override = default;

private slots:
    void digestProteinTest();
    void digestProteinTestTestPartial();
    void ddigestProteinTestTestRagged();
    void proteomicsTest();
};


void ProteinDigestomaticTests::digestProteinTest()
{
    ERR_INIT;

    const QString proteinSequence = QStringLiteral("QNITEEFYQSTCSAVSKASDF");
    const QString proteinSequence2 = QStringLiteral("QNITEEFYQSTCSAVSKASDFK");

    PythiaParameters params;
    params.cTermCleavePoints << "K";
    params.nTermCleavePoints << "Q";
    params.allowedMissedCleavages = 0;
    params.raggedness = Raggedness::NoRagged;
    params.peptideLengthMin = 1;
    params.peptideLengthMax = 1000;

    ProteinDigestomatic proteinDigestomatic(params);

    QVector<PeptideSequence> peptideSequences;
    e = proteinDigestomatic.digestProtein(
            proteinSequence,
            &peptideSequences
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideSequences.size(), 3);

    const PeptideSequence &first = peptideSequences.front();
    const PeptideSequence &mid = peptideSequences.at(1);
    const PeptideSequence &last = peptideSequences.back();

    QCOMPARE(first.startIndex, 0);
    QCOMPARE(first.endIndex, 7);
    QCOMPARE(first.previousResidue, '*');
    QCOMPARE(first.firstResidue, 'Q');
    QCOMPARE(first.lastResidue, 'Y');
    QCOMPARE(first.postResidue, 'Q');
    QCOMPARE(first.sequence, QStringLiteral("QNITEEFY"));

    QCOMPARE(mid.startIndex, 8);
    QCOMPARE(mid.endIndex, 16);
    QCOMPARE(mid.previousResidue, 'Y');
    QCOMPARE(mid.firstResidue, 'Q');
    QCOMPARE(mid.lastResidue, 'K');
    QCOMPARE(mid.postResidue, 'A');
    QCOMPARE(mid.sequence, QStringLiteral("QSTCSAVSK"));

    QCOMPARE(last.startIndex, 17);
    QCOMPARE(last.endIndex, 20);
    QCOMPARE(last.previousResidue, 'K');
    QCOMPARE(last.firstResidue, 'A');
    QCOMPARE(last.lastResidue, 'F');
    QCOMPARE(last.postResidue, '*');
    QCOMPARE(last.sequence, QStringLiteral("ASDF"));

    peptideSequences.clear();
    e = proteinDigestomatic.digestProtein(
            proteinSequence2,
            &peptideSequences
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(peptideSequences.size(), 3);

    const PeptideSequence &last2 = peptideSequences.back();

    QCOMPARE(last2.startIndex, 17);
    QCOMPARE(last2.endIndex, 21);
    QCOMPARE(last2.previousResidue, 'K');
    QCOMPARE(last2.firstResidue, 'A');
    QCOMPARE(last2.lastResidue, 'K');
    QCOMPARE(last2.postResidue, '*');
    QCOMPARE(last.sequence, QStringLiteral("ASDFK"));

    e = proteinDigestomatic.digestProtein(
            QString(),
            &peptideSequences
            );

    QCOMPARE(e, eError);
}

void ProteinDigestomaticTests::digestProteinTestTestPartial()
{
    ERR_INIT;

    const QString proteinSequence = QStringLiteral("AAKCCKDDKEE");

    PythiaParameters params;
    params.cTermCleavePoints << "K";
    params.allowedMissedCleavages = 2;
    params.raggedness = Raggedness::NoRagged;
    params.peptideLengthMin = 1;
    params.peptideLengthMax = 1000;


    ProteinDigestomatic proteinDigestomatic(params);

    QVector<PeptideSequence> peptideSequences;
    e = proteinDigestomatic.digestProtein(proteinSequence, &peptideSequences);
    QCOMPARE(e, eNoError);

    const QStringList expected = {
            "AAK",
            "CCK",
            "DDK",
             "EE",
            "AAKCCK",
            "CCKDDK",
             "DDKEE",
            "AAKCCKDDK",
             "CCKDDKEE"
    };

    QCOMPARE(peptideSequences.size(), expected.size());

    for(int i = 0; i < expected.size(); i++){
        QVERIFY(expected.contains(peptideSequences.at(i).sequence));
    }

}

void ProteinDigestomaticTests::ddigestProteinTestTestRagged()
{
    ERR_INIT;

    const QString proteinSequence = QStringLiteral("ACDKEFGKHILKMNP");

    PythiaParameters params;
    params.cTermCleavePoints << "K";
    params.allowedMissedCleavages = 0;
    params.raggedness = Raggedness::BothRagged;
    params.peptideLengthMin = 1;
    params.peptideLengthMax = 1000;

    ProteinDigestomatic proteinDigestomatic(params);

    QVector<PeptideSequence> peptideSequences;
    e = proteinDigestomatic.digestProtein(proteinSequence, &peptideSequences);

    QCOMPARE(e, eNoError);

    const QStringList expected = {
             "ACDK",
             "EFGK",
             "HILK",
             "MNP",
             "ACD",
             "CDK",
             "AC",
             "DK",
             "EFG",
             "FGK",
             "EF",
             "GK",
             "HIL",
             "ILK",
             "HI",
             "LK",
             "MN",
             "NP"
    };



    for(int i = 0; i < expected.size(); i++){
        QCOMPARE(peptideSequences.at(i).sequence, expected.at(i));
    }

    const PeptideSequence &details = peptideSequences.at(7);
    QCOMPARE(details.startIndex, 2);
    QCOMPARE(details.endIndex, 3);
    QCOMPARE(details.previousResidue, 'C');
    QCOMPARE(details.firstResidue, 'D');
    QCOMPARE(details.lastResidue, 'K');
    QCOMPARE(details.postResidue, 'E');
}

void ProteinDigestomaticTests::proteomicsTest()
{
    //TODO Change hard coded path to use qapp->directory.
    const QString &fastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    QElapsedTimer et;
    et.start();
    FastaReader fastaReader;
    fastaReader.parseFastaFile(fastaFilePath);
    qDebug() << QStringLiteral("Fasta read in:") << et.elapsed() / 1000.0;
    
    ERR_INIT
    // TODO correct this.
    PythiaParameters params;
    params.cTermCleavePoints << "K" << "R" << "F" << "L"  << "W" ;
    params.nTermCleavePoints << "D" << "E";
    params.allowedMissedCleavages = 1;
    params.raggedness = Raggedness::NTermRagged;
    params.peptideLengthMin = 6;
    params.peptideLengthMax = 40;

    ProteinDigestomatic proteinDigestomatic(params);
    const QVector<FastaEntry> proteins = fastaReader.fastaEntries().values().toVector();

    QElapsedTimer et2;
    et2.start();
    int counter = 0;
    int peptideCount = 0;
    for(const FastaEntry &fe : proteins){

        QVector<PeptideSequence> peptideSequences;
        e = proteinDigestomatic.digestProtein(fe.fastaSequence, &peptideSequences);

        peptideCount += peptideSequences.size();
        if(++counter % 10000 == 0)
            qDebug() << counter;
    }

    const double timeElapsedSeconds = et2.elapsed() / 1000.0;
    qDebug() << peptideCount << timeElapsedSeconds;
    QCOMPARE(peptideCount, 338191); //TODO updated test
    QVERIFY(timeElapsedSeconds < 20);
}

QTEST_MAIN(ProteinDigestomaticTests)

#include "ProteinDigestomaticTests.moc"
