#include "LibraryBuilderWorkFlow.h"

#include "PythiaParameterReader.h"
#include "Error.h"

#include <QtTest>

class LibraryBuilderWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    LibraryBuilderWorkFlowTests() = default;
    ~LibraryBuilderWorkFlowTests() override = default;

private slots:

    void execTest();


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
            fastaFilePath()
            );

    QCOMPARE(e, eNoError);

}


QTEST_MAIN(LibraryBuilderWorkFlowTests)

#include "LibraryBuilderWorkFlowTests.moc"
