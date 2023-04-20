#include "LibraryBuilderWorkFlow.h"

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

        const QString &paramsFile
                = QDir(qApp->applicationDirPath()).filePath("WorkFlowTestsParams.pythia");

        PythiaParameterReader reader;
        PythiaParameters pythiaParameters;
        reader.readFile(paramsFile);
        reader.loadPythiaParameters(&pythiaParameters);

        pythiaParameters.topNMs2Ions = 12;
//        pythiaParameters.ms2ExtractionWidthPPM = 20;

        return pythiaParameters;
    }

};

void LibraryBuilderWorkFlowTests::execTest() {

    ERR_INIT

    const QString peptidesCSVFilePath = QDir(qApp->applicationDirPath()).filePath("peptides.csv");
    
    const QString model1FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_1.hdf5.json");
    const QString model2FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_2.hdf5.json");
    const QString model3FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_3.hdf5.json");
    const QString model4FilePath
            = QDir(qApp->applicationDirPath()).filePath("rnn_linear_charge_w_precursors_nce_4.hdf5.json");

    QString fragLibFilePath;

    LibraryBuilderWorkFlow libraryBuilderWorkFlow;
    e = libraryBuilderWorkFlow.exec(
            peptidesCSVFilePath,
            &fragLibFilePath
            );
    QCOMPARE(e, eError);

    e = libraryBuilderWorkFlow.init(
            pythiaParameters(),
            model1FilePath,
            model2FilePath,
            model3FilePath,
            model4FilePath
            );
    QCOMPARE(e, eNoError);

    e = libraryBuilderWorkFlow.exec(
            peptidesCSVFilePath,
            &fragLibFilePath
            );
    QCOMPARE(e, eNoError);

    QFileInfo fi(fragLibFilePath);
    QCOMPARE(fi.exists(), true);

}


QTEST_MAIN(LibraryBuilderWorkFlowTests)

#include "LibraryBuilderWorkFlowTests.moc"
