#include "SequenceSubstringSearchomatic.h"

#include "Error.h"
#include "FastaReader.h"

#include <QtTest>

using namespace Error;

class SequenceSubstringSearchomaticTests : public QObject
{
    Q_OBJECT

public:
    SequenceSubstringSearchomaticTests() = default;
    ~SequenceSubstringSearchomaticTests() override = default;

private slots:

    static void searchSequencesTest();

};

void SequenceSubstringSearchomaticTests::searchSequencesTest() {
    ERR_INIT

    const QString &fastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath);
    QCOMPARE(e, eNoError);

    SequenceSubstringSearchomatic sequenceSubstringSearchomatic;
    e = sequenceSubstringSearchomatic.init(fastaReader.fastaEntries());
    QCOMPARE(e, eNoError);

    const QString searchSequence = QStringLiteral("FLEEDGNVNSKLTKDSVWNYHC");

    QVector<QString> proteinGroups;
    e = sequenceSubstringSearchomatic.findProteinGroups(
        {searchSequence},
        &proteinGroups
        );
    QCOMPARE(e, eNoError);

    const QString expectedResult
        = QStringLiteral(">sp|P00488|F13A_HUMAN Coagulation factor XIII A chain OS=Homo sapiens OX=9606 GN=F13A1 PE=1 SV=4");

    QCOMPARE(proteinGroups.front(), expectedResult);

}


QTEST_MAIN(SequenceSubstringSearchomaticTests)

#include "SequenceSubstringSearchomaticTests.moc"
