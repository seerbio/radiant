//
// Created by anichols on 7/16/24.
//

#include "AWSStreamifier.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <iostream>
#include <sstream>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <tuple>
#include <string>

#include <QDebug>

// #include "ErrorUtils.h"


bool AWSStreamifier::setAWSCredentials(
    const QString& accessKeyId,
    const QString& secretAccessKey,
    const QString& accessToken
    ) {

    m_accessKeyId = accessKeyId;
    m_secretAccessKey = secretAccessKey;
    m_accessToken = accessToken;

    if (!credentialsValid()) {
        return false;
    }



    return true;
}

bool AWSStreamifier::credentialsValid() const {
    return !m_accessToken.isEmpty() && !m_secretAccessKey.isEmpty() && !m_accessToken.isEmpty();
}

namespace {

    std::tuple<bool, QString, QString> getBucketNameVsObjectKey(const QString &uri) {

        const QString s3Hash = QStringLiteral("s3://");

        QString uriCopy = uri;
        uriCopy = uriCopy.replace("\\", "/");
        uriCopy = uriCopy.replace(s3Hash, "");

        QStringList uriCopySplit = uriCopy.split("/");

        constexpr int minSplitSize = 2;
        if (uriCopySplit.size() < minSplitSize) {
            return {false, {}, {}};
        }

        const QString bucketName = uriCopySplit.at(0);

        uriCopySplit.pop_front();
        const QString objectkey = uriCopySplit.join("/");

        return {true, bucketName, objectkey};
    }

    QPair<bool, Aws::IOStream&> getAWSStream(
        const QString &accessKeyId,
        const QString &secretAccessKey,
        const QString &accessToken,
        const Aws::String &bucketName,
        const Aws::String &objectKey
    ) {

        Aws::SDKOptions options;
        Aws::InitAPI(options);
        {
            Aws::Client::ClientConfiguration clientConfiguration;
            clientConfiguration.region = Aws::Region::US_WEST_2;

            const Aws::Auth::AWSCredentials credentials(
                accessKeyId.toStdString(),
                secretAccessKey.toStdString(),
                accessToken.toStdString()
                );

            Aws::S3::S3Client s3Client(
                credentials,
                nullptr,
                clientConfiguration
                );

            Aws::S3::Model::GetObjectRequest objectRequest;
            objectRequest.WithBucket(bucketName).WithKey(objectKey);

            const auto objectOutcome = s3Client.GetObject(objectRequest);

            if (objectOutcome.IsSuccess()) {

                Aws::IOStream& s3Stream = objectOutcome.GetResult().GetBody();
                return {true, s3Stream};

            } else {
                qCritical() << "Failed to get object:";
                qCritical() << QString::fromStdString(objectOutcome.GetError().GetMessage());
            }
        }
        Aws::ShutdownAPI(options);
    }

}//namespace
bool AWSStreamifier::streamParquetFile(
    const QString& uri,
    std::vector<uint8_t>* fileBuffer
    ) {

    if (!credentialsValid()) {
        return false;
    }

    return true;
}

bool AWSStreamifier::streamTextFile(
    const QString& uri,
    char** fileData
    ) {

    if (!credentialsValid()) {
        return false;
    }

    const std::tuple<bool, QString, QString> bucketNameVsObjectKey = getBucketNameVsObjectKey(uri);

    const bool isValid = std::get<0>(bucketNameVsObjectKey);
    if (!isValid) {
        return false;
    }

    const QString bucket = std::get<1>(bucketNameVsObjectKey);
    const QString objectKey = std::get<2>(bucketNameVsObjectKey);;

    Aws::String accessKeyId = "ASIA2OXTCKO6EKHAAK65";
    Aws::String secretAccessKey = "hsMVg/bbS1XiTz8+hVcaGV95XOGrJ3+vqrNn27w9";
    Aws::String accessToken = "IQoJb3JpZ2luX2VjEDQaCXVzLXdlc3QtMiJHMEUCIB1X6orkv5E/ho3KWXWFV59MxmkgUngISc3gF+yvzc+sAiEAogdT5AlcQqLPlrNvKDBLeiNt4g+rbEkTMbtNm9/cVXYqlQMI/f//////////ARAFGgw3MTg4NDMwNDA3MDAiDA5/vIUPXpCAvkRrGirpAhIaGD4capupKRYXoKJz9T3/xsHyEz1BUf1TFrCo2LtGMpoFfnjr0/fiQMs6M74ZFSa+F84ImhZRpHYMOtSIWrslFbquWB4OGZhdvSOcVU3Dn3qTccTXIwlYbCVCT30KSdbLFmzzHvf9iIhOhwU4tSvLXvGQEZTWAlI+8h34oQ7uKXfYrhFnxe2NbQFTjJDD2Y5TtDC+mTN8wYRIHvR2HcYqtcgxIr+kJfZ0kWyzPfNcbpLTTiFu3GB7F0Olr89MY/XMn6b6zuH5YSgu5vLVuD9/oGl2CT9W/fVEqgrxpWRgd6tOy00JGZALhX8dCmUqePXGkfpUeUHdGY1QIN9U8Rb64nxksN33/VbKrLt9H/dkiLLt4rJz7tCWXQqgTUOp3QaA3nlJPRrTX+uxKqXIGtajvsFbEno/k/OtifK/up7KewO7Q/+7TpLBBQn5KpqA98otZ59pI92eXYK6BmACG9Ygd8EvVtcey/IwnqPitAY6pgHwLgoxih2LYYH5J9yyKGaYDL8Z6UZQV7U/yuLe1oUV8wU3ln/+zlPXNTyXvkw//m297PwZHH5LbqI0355Vu9eN2GmzCkU1YdyxUHezKbyp3SdvM0DF0+yNH1GHoWne4Qqoy2uoOcPSqh0HUaS/RI+u9K5jrXEg9dl085XBM3DItpwGRZCv2lg887zmX91/iJwYzfuFuQsjybVYfzerqXRbNBexgNn7";

    Aws::SDKOptions options;
    Aws::InitAPI(options);
    {
        Aws::Client::ClientConfiguration clientConfiguration;
        clientConfiguration.region = Aws::Region::US_WEST_2;

        Aws::Auth::AWSCredentials credentials(
            m_accessKeyId.toStdString(),
            m_secretAccessKey.toStdString(),
            m_accessToken.toStdString()
            );

        Aws::S3::S3Client s3_client(credentials, nullptr, clientConfiguration);

        Aws::String bucket_name = "seer-experiments";
        Aws::String object_key = "EXP24007/mzml/EXP24007_2024ms0113X41_A.raw.mzML";
        object_key = "EXP24001/mzml/EXP24001_2023ms1005bX1_A.raw.mzML";

        Aws::S3::Model::GetObjectRequest object_request;
        object_request.WithBucket(bucket_name).WithKey(object_key);

        auto get_object_outcome = s3_client.GetObject(object_request);

        if (get_object_outcome.IsSuccess()) {

            Aws::IOStream& s3_stream = get_object_outcome.GetResult().GetBody();

            std::string line;
            while (std::getline(s3_stream, line)) {
                std::cout << line << std::endl;
            }

            // Now you can use `buffer` as needed
            s3_stream.clear();
            s3_stream.flush();

        } else {
            // Handle errors
            std::cerr << "Failed to get object: " <<
                      get_object_outcome.GetError().GetMessage() << std::endl;
        }
    }
    Aws::ShutdownAPI(options);

    // QPair<bool, Aws::IOStream&> ioStreamResult = getAWSStream(
    //     m_accessKeyId,
    //     m_secretAccessKey,
    //     m_accessToken,
    //     bucket.toStdString(),
    //     objectKey.toStdString()
    //     );
    //
    // // if (!ioStreamResult.first) {
    // //     return false;
    // // }


    return true;
}
