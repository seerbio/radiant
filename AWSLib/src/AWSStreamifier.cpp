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

}//namespace
bool AWSStreamifier::streamParquetFile(
    const QString& uri,
    std::vector<uint8_t>* fileBuffer
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
    const QString objectKey = std::get<2>(bucketNameVsObjectKey);

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

        Aws::S3::S3Client s3Client(credentials, nullptr, clientConfiguration);

        Aws::S3::Model::GetObjectRequest objectRequest;
        objectRequest.WithBucket(bucket.toStdString()).WithKey(objectKey.toStdString());

        auto getObjectOutcome = s3Client.GetObject(objectRequest);

        if (getObjectOutcome.IsSuccess()) {

            auto& s3Stream = getObjectOutcome.GetResult().GetBody();
            std::vector<uint8_t> buffer(std::istreambuf_iterator<char>(s3Stream), {});
            fileBuffer = &buffer;

            std::cout << "Successfully read Parquet data into memory." << std::endl;
        } else {
            std::cerr << "Failed to get object: " << getObjectOutcome.GetError().GetMessage() << std::endl;
        }

    }
    Aws::ShutdownAPI(options);

    return true;
}

bool AWSStreamifier::streamTextFile(
    const QString& uri,
    uchar** fileData
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
    const QString objectKey = std::get<2>(bucketNameVsObjectKey);

    std::string lineAggregate;

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

        Aws::S3::S3Client s3Client(credentials, nullptr, clientConfiguration);

        Aws::S3::Model::GetObjectRequest objectRequest;
        objectRequest.WithBucket(bucket.toStdString()).WithKey(objectKey.toStdString());

        auto getObjectOutcome = s3Client.GetObject(objectRequest);

        if (getObjectOutcome.IsSuccess()) {

            Aws::IOStream& s3Stream = getObjectOutcome.GetResult().GetBody();

            std::string line;
            while (std::getline(s3Stream, line)) {
                lineAggregate += line;
            }

            s3Stream.clear();
            s3Stream.flush();

        } else {
            std::cerr << "Failed to get object: " <<
                      getObjectOutcome.GetError().GetMessage() << std::endl;
        }
    }
    Aws::ShutdownAPI(options);

    auto dataPtr = reinterpret_cast<uchar*>(lineAggregate.data());
    fileData = &dataPtr;

    return true;
}
