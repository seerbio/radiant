//
// Created by Drucifer on 12/4/2021.
//

#include "JsonParametersReader.h"

#include "ErrorUtils.h"

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>

bool hasJsonError(const QJsonParseError &parseError) {

    if (parseError.error != QJsonParseError::NoError) {
        qDebug() << "Parsing error:" << parseError.errorString();
        return true;
    }

    return false;
}


Err JsonParametersReader::readFile(const QString &filePath) {

    ERR_INIT

    QFile file(filePath);
    const bool fileOpened = file.open(QIODevice::ReadOnly | QIODevice::Text);
    e = ErrorUtils::isEqual(fileOpened, true, Error::eFileError); ree;

    const QString jsonString = file.readAll();
    file.close();

    QJsonParseError parseError;
    QJsonDocument document = QJsonDocument::fromJson(jsonString.toUtf8(), &parseError);
    e = ErrorUtils::isFalse(hasJsonError(parseError)); ree;


    QJsonObject jsonObject = document.object();

    for(const QString& key : jsonObject.keys()){
        const QVariant value = jsonObject.value(key).toVariant();
        m_jsonContents.insert(key, value);
    }

    e = ErrorUtils::isNotEmpty(m_jsonContents); ree;

    ERR_RETURN
}

QMap<QString, QVariant> JsonParametersReader::JsonParametersReader::jsonContents() const {
    return m_jsonContents;
}

bool JsonParametersReader::jsonContentsIsEmpty() const {
    return m_jsonContents.isEmpty();
}
