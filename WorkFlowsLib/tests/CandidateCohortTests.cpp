#include "PythiaDIAFFWorkflow.h"
#include "ParquetReader.h"

#include <QFile>
#include <QTemporaryDir>
#include <QtTest>

class CandidateCohortTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void preparedLibraryIsBoundToItsSourceAndDecoyMode();
    void rejectsInvalidViewsBeforeOpeningInput();
};

void CandidateCohortTests::preparedLibraryIsBoundToItsSourceAndDecoyMode() {
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    FragLibReaderRow row;
    row.peptideSequenceChargeKey = "ACDEFGHIK|2";
    row.precursorCharge = 2;
    row.iRT = 1.0f;
    row.mass = 1000.0;
    row.mzVals = {300.0f, 400.0f, 500.0f, 600.0f};
    row.intensityVals = {1.0f, .8f, .6f, .4f};
    row.ionLabels = "b2;b3;y2;y3";
    row.proteinGroups = "P1";
    const QString libraryPath = temporary.filePath("test.fragLibFF");
    QCOMPARE(ParquetReader::write(QList<FragLibReaderRow>{row}, libraryPath), eNoError);
    PythiaDIAFFWorkflow::LibraryHandle library;
    QCOMPARE(PythiaDIAFFWorkflow::prepareLibrary(libraryPath, false, &library), eNoError);
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.threadCount = 1;
    parameters.peptideLengthMin = 7;
    parameters.useAlternativeDecoys = false;
    PythiaDIAFFWorkflow first, second;
    QCOMPARE(first.init(parameters, libraryPath, {}, {}, library), eNoError);
    QCOMPARE(second.init(parameters, libraryPath, {}, {}, library), eNoError);
    auto changed = parameters;
    changed.useAlternativeDecoys = true;
    PythiaDIAFFWorkflow wrongMode;
    QCOMPARE(wrongMode.init(changed, libraryPath, {}, {}, library), eValueError);
    const QString otherPath = temporary.filePath("other.fragLibFF");
    QVERIFY(QFile::copy(libraryPath, otherPath));
    PythiaDIAFFWorkflow wrongSource;
    QCOMPARE(wrongSource.init(parameters, otherPath, {}, {}, library), eValueError);
    QFile modified(libraryPath);
    QVERIFY(modified.open(QIODevice::Append));
    QCOMPARE(modified.write("x", 1), qint64(1));
    modified.close();
    PythiaDIAFFWorkflow changedFile;
    QCOMPARE(changedFile.init(parameters, libraryPath, {}, {}, library), eValueError);
}

void CandidateCohortTests::rejectsInvalidViewsBeforeOpeningInput() {
    PythiaDIAFFWorkflow workflow;
    // Ordinary report mode must not silently become a candidate-only run.
    QCOMPARE(workflow.processCandidateViews("missing.mzML", {{4, 4, "result"}}), eValueError);
    QCOMPARE(workflow.processCandidateViews("missing.mzML", {}), eValueError);
}

QTEST_MAIN(CandidateCohortTests)
#include "CandidateCohortTests.moc"
