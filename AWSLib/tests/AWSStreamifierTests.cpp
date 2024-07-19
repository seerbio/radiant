//
// Created by anichols on 11/07/2021.
//

#include "AWSStreamifier.h"

#include <QProcessEnvironment>
#include <QtTest/QtTest>

class AWSStreamifierTests : public QObject
{
    Q_OBJECT

public:
    AWSStreamifierTests() = default;
    ~AWSStreamifierTests() override = default;

private Q_SLOTS:

    void streamTextFileTest();

    void streamParquetFileTest();
};

void AWSStreamifierTests::streamTextFileTest() {

    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();

    if(
        !(env.contains("AWS_ACCESS_KEY_ID")
        && env.contains("AWS_SECRET_ACCESS_KEY")
        && env.contains("AWS_SESSION_TOKEN"))
        ) {
        QSKIP("Skipping due to no credentials provided");
    }

    AWSStreamifier awsStreamifier;
    const bool credentialsSet = awsStreamifier.setAWSCredentials(
        env.value("AWS_ACCESS_KEY_ID"),
        env.value("AWS_SECRET_ACCESS_KEY"),
        env.value("AWS_SESSION_TOKEN")
        );

    const QString uri = "s3://seer-experiments/EXP24001/mzml/EXP24001_2023ms1005bX1_A.raw.mzML";

    uchar* fileData;
    const bool streamIsValid = awsStreamifier.streamTextFile(uri, &fileData);

}

void AWSStreamifierTests::streamParquetFileTest() {

    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();

    if(
        !(env.contains("AWS_ACCESS_KEY_ID")
        && env.contains("AWS_SECRET_ACCESS_KEY")
        && env.contains("AWS_SESSION_TOKEN"))
        ) {
        QSKIP("Skipping due to no credentials provided");
        }

    AWSStreamifier awsStreamifier;
    const bool credentialsSet = awsStreamifier.setAWSCredentials(
        env.value("AWS_ACCESS_KEY_ID"),
        env.value("AWS_SECRET_ACCESS_KEY"),
        env.value("AWS_SESSION_TOKEN")
        );

    const QString uri = "s3://seer-experiments/EXP24039/interimresult/0b6d21d0b237173e9ed0ec96aea642cd/500873.EXP24039_2024ms0407dX45_A.raw.mzML.pythiaDIA";

    std::vector<uint8_t> fileBuffer;
    const bool streamIsValid = awsStreamifier.streamParquetFile(uri, &fileBuffer);

}

QTEST_MAIN(AWSStreamifierTests)
#include "AWSStreamifierTests.moc"
