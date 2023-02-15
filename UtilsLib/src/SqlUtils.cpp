//
// Created by Drucifer on 11/21/2021.
//

#include "SqlUtils.h"

#include <QCryptographicHash>
#include <QFile>


QByteArray SqlUtils::fileChecksum(const QString &filePath) {

    QFile fileName(filePath);

    if (fileName.open(QFile::ReadOnly)) {
        QCryptographicHash hash(QCryptographicHash::Md5);
        if (hash.addData(&fileName)) {
            return hash.result();
        }
    }
    return {};
}
