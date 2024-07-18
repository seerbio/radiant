//
// Created by anichols on 7/16/24.
//

#ifndef AWSSTREAMIFIER_H
#define AWSSTREAMIFIER_H


#include "AWSLib_Exports.h"

/* NOTE:
 * Cannot use Error.h or ErrorUtils.h as includes as it conflict w/ Arrow/Parquet libary.
 */

#include <QString>

class AWSLIB_EXPORTS AWSStreamifier {

public:

    AWSStreamifier() = default;
    ~AWSStreamifier() = default;

    bool setAWSCredentials(
        const QString &accessKeyId,
        const QString &secretAccessKey,
        const QString &accessToken
        );

    bool streamParquetFile(
        const QString &uri,
        std::vector<uint8_t> *fileBuffer
        );

    bool streamTextFile(
        const QString &uri,
        char* *fileData
        );

private:

    [[nodiscard]] bool credentialsValid() const;

private:

    QString m_accessKeyId;
    QString m_secretAccessKey;
    QString m_accessToken;

};



#endif //AWSSTREAMIFIER_H
