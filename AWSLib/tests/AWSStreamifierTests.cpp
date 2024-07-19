//
// Created by anichols on 11/07/2021.
//

#include "AWSStreamifier.h"

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

    QSKIP("temp");

    const QString accessKeyId = "ASIA2OXTCKO6AFXJIBOJ";
    const QString secretAccessKey = "nyAGcu1/gHDBqypr2UZQVl1RNTa5VA8joOz4hC+E";
    const QString accessToken = "IQoJb3JpZ2luX2VjEFsaCXVzLXdlc3QtMiJGMEQCICxWGVgGZSdqDtcHBpIh5zyHB3xNjOvNjQ5w0t8W+n8pAiAGt/oyjpFGa2LFjKI8vsNCpngZaFpplDl5aHWkVLMc9SqMAwg0EAUaDDcxODg0MzA0MDcwMCIMjQq949PgLVNtnPTKKukCES+e7JHJST1dkSowNU0HbRbe0jtEk8lViyUNqTvEhvo5Z/Bs7RDLgX7XoWT6Wnmcx2creKM9PTMUP8Yi7HrIC6jxCG1DQeJpxllwFd3dT1KEk0rs4VKrFkDXeP3cmwIHhGuoQC9EM6tDfYG2uMtskVIdIOQ18DJvLdI4UksqIC4iZu1ADF/6ctOmVUnOeSmhM9eGHtHd3mrPiw15wSG/gQFYQ5rUy+vbO4yCm1kzSUP9aVxkZpl6K6uMlv1d38fGhbe9OIZncj7Ozt1CuquP9XJBV1DnTEl9YDmNjii+3Ck5+t429tJ3rES3QHBxr4a1Somt6ZoSysgnyogr5qbNrfD80EmEa/werYG9edIy6fGc+r5SAe10k/AfzOO6nHtSDKVzUquRwEmHxls/vXAWz1LTjJuQ5Zap1YcSW2Ob7w+u6H0rNwRC1xgNJIT8+4kLgywKToiCbFVh9R1OXp8C3LyBDCE+CoYH/DD45+q0BjqnAc8Qa9F7FlWe7NR6Wz1+/RzH11O5gBQ7mLYd+oW3WhEwR0TRndRRD9X812sSK43/3gQG7sEDNXrHYRIcCzoDcu3wayGa8KXTiLbMHIXwrtmBa3kiPX54HvNGCPPceMYhWoozejRMgxMWsLSN8RclngQDy9Fhy6gvCm5r2qYghZk0L02xBRpECpfo+x8ZGZ5/adpuAh2etw8bVTmeVlN+cTJFwOHWaA4S";

    AWSStreamifier awsStreamifier;
    const bool credentialsSet = awsStreamifier.setAWSCredentials(
        accessKeyId,
        secretAccessKey,
        accessToken
        );

    const QString uri = "s3://seer-experiments/EXP24001/mzml/EXP24001_2023ms1005bX1_A.raw.mzML";

    uchar* fileData;
    const bool streamIsValid = awsStreamifier.streamTextFile(uri, &fileData);

}

void AWSStreamifierTests::streamParquetFileTest() {

    // QSKIP("temp");

    const QString accessKeyId = "ASIA2OXTCKO6AFXJIBOJ";
    const QString secretAccessKey = "nyAGcu1/gHDBqypr2UZQVl1RNTa5VA8joOz4hC+E";
    const QString accessToken = "IQoJb3JpZ2luX2VjEFsaCXVzLXdlc3QtMiJGMEQCICxWGVgGZSdqDtcHBpIh5zyHB3xNjOvNjQ5w0t8W+n8pAiAGt/oyjpFGa2LFjKI8vsNCpngZaFpplDl5aHWkVLMc9SqMAwg0EAUaDDcxODg0MzA0MDcwMCIMjQq949PgLVNtnPTKKukCES+e7JHJST1dkSowNU0HbRbe0jtEk8lViyUNqTvEhvo5Z/Bs7RDLgX7XoWT6Wnmcx2creKM9PTMUP8Yi7HrIC6jxCG1DQeJpxllwFd3dT1KEk0rs4VKrFkDXeP3cmwIHhGuoQC9EM6tDfYG2uMtskVIdIOQ18DJvLdI4UksqIC4iZu1ADF/6ctOmVUnOeSmhM9eGHtHd3mrPiw15wSG/gQFYQ5rUy+vbO4yCm1kzSUP9aVxkZpl6K6uMlv1d38fGhbe9OIZncj7Ozt1CuquP9XJBV1DnTEl9YDmNjii+3Ck5+t429tJ3rES3QHBxr4a1Somt6ZoSysgnyogr5qbNrfD80EmEa/werYG9edIy6fGc+r5SAe10k/AfzOO6nHtSDKVzUquRwEmHxls/vXAWz1LTjJuQ5Zap1YcSW2Ob7w+u6H0rNwRC1xgNJIT8+4kLgywKToiCbFVh9R1OXp8C3LyBDCE+CoYH/DD45+q0BjqnAc8Qa9F7FlWe7NR6Wz1+/RzH11O5gBQ7mLYd+oW3WhEwR0TRndRRD9X812sSK43/3gQG7sEDNXrHYRIcCzoDcu3wayGa8KXTiLbMHIXwrtmBa3kiPX54HvNGCPPceMYhWoozejRMgxMWsLSN8RclngQDy9Fhy6gvCm5r2qYghZk0L02xBRpECpfo+x8ZGZ5/adpuAh2etw8bVTmeVlN+cTJFwOHWaA4S";

    AWSStreamifier awsStreamifier;
    const bool credentialsSet = awsStreamifier.setAWSCredentials(
        accessKeyId,
        secretAccessKey,
        accessToken
        );

    const QString uri = "s3://seer-experiments/EXP24039/interimresult/0b6d21d0b237173e9ed0ec96aea642cd/500873.EXP24039_2024ms0407dX45_A.raw.mzML.pythiaDIA";

    std::vector<uint8_t> fileBuffer;
    const bool streamIsValid = awsStreamifier.streamParquetFile(uri, &fileBuffer);

}

QTEST_MAIN(AWSStreamifierTests)
#include "AWSStreamifierTests.moc"
