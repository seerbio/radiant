//
// Created by Drucifer on 12/4/2021.
//

#ifndef JSONPARAMETERSREADER_H
#define JSONPARAMETERSREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"

#include <QString>
#include <QVariant>

using namespace Error;

class FILEREADERSLIB_EXPORTS JsonParametersReader {

public:
    JsonParametersReader() = default;
    ~JsonParametersReader() = default;

    Err readFile(const QString &filePath);

    [[nodiscard]] QMap<QString, QVariant> jsonContents() const;


protected:

    [[nodiscard]] bool jsonContentsIsEmpty() const;


private:

    QMap<QString, QVariant> m_jsonContents;



};


#endif //JSONPARAMETERSREADER_H
