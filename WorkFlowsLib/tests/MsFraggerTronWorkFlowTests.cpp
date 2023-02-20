#include "MsFraggerTronWorkFlow.h"

#include "Error.h"
#include "PeptidesLibraryTron.h"
#include "PythiaParameterReader.h"

#include <QtTest>

class MsFraggerTronWorkFlowTests : public QObject
{
    Q_OBJECT

public:

    MsFraggerTronWorkFlowTests() = default;
    ~MsFraggerTronWorkFlowTests() override = default;

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

};

void MsFraggerTronWorkFlowTests::execTest() {

    ERR_INIT

    MsFraggerTronWorkFlow msFraggerTronWorkFlow(pythiaParameters());

    const QString mzMLFileURI
        = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

    e = msFraggerTronWorkFlow.exec(mzMLFileURI);

    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsFraggerTronWorkFlowTests)

#include "MsFraggerTronWorkFlowTests.moc"
