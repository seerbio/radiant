//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "FastaReader.h"

#include <QtTest/QtTest>

class FastaReaderTests : public QObject
{
    Q_OBJECT

public:
    FastaReaderTests() = default;
    ~FastaReaderTests() override = default;

private Q_SLOTS:

    void parseFastaFileTest();
    void hasDecoysTest();


};

void FastaReaderTests::parseFastaFileTest() {

    ERR_INIT;

    const QString &fastaFilePath
        = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    FastaReader fastaReader;
    fastaReader.parseFastaFile(fastaFilePath, false);

    const QVector<FastaEntry> proteins = fastaReader.fastaEntries().values().toVector();

    const QString expectedFastaDescription
        = QStringLiteral(">sp|O60884|DNJA2_HUMAN DnaJ homolog subfamily A member 2 OS=Homo sapiens OX=9606 GN=DNAJA2 PE=1 SV=1");

    const QString expectedFastaSequence
        = QStringLiteral("MANVADTKLYDILGVPPGASENELKKAYRKLAKEYHPDKNPNAGDKFK"
                         "EISFAYEVLSNPEKRELYDRYGEQGLREGSGGGGGMDDIFSHIFGGGL"
                         "FGFMGNQSRSRNGRRRGEDMMHPLKVSLEDLYNGKTTKLQLSKNVLCSA"
                         "CSGQGGKSGAVQKCSACRGRGVRIMIRQLAPGMVQQMQSVCSDCNGEGE"
                         "VINEKDRCKKCEGKKVIKEVKILEVHVDKGMKHGQRITFTGEADQAPGVEP"
                         "GDIVLLLQEKEHEVFQRDGNDLHMTYKIGLVEALCGFQFTFKHLDGRQIVVK"
                         "YPPGKVIEPGCVRVVRGEGMPQYRNPFEKGDLYIKFDVQFPENNWINPDKL"
                         "SELEDLLPSRPEVPNIIGETEEVELQEFDSTRGSGGGQRREAYNDSSDEESSSH"
                         "HGPGVQCAHQ");

    QCOMPARE(proteins.size(), 996);
    QCOMPARE(proteins.at(0).fastaDescription, expectedFastaDescription);
    QCOMPARE(proteins.at(0).fastaSequence, expectedFastaSequence);



    FastaReader fastaReaderDecoys;
    fastaReaderDecoys.parseFastaFile(fastaFilePath, true);

    const QVector<FastaEntry> proteinsDecoys = fastaReaderDecoys.fastaEntries().values().toVector();
    QCOMPARE(proteinsDecoys.size(), 1992);

}

void FastaReaderTests::hasDecoysTest() {

    //TODO Change hard coded path to use qapp->directory.
    const QString &fastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    FastaReader fastaReader;
    fastaReader.parseFastaFile(fastaFilePath, false);

    QCOMPARE(fastaReader.hasDecoys(), false);

}

QTEST_MAIN(FastaReaderTests)
#include "FastaReaderTests.moc"
