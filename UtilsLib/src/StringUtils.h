#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include "UtilsLib_Exports.h"

#include <QDebug>
#include <QMap>
#include <QMetaEnum>
#include <QString>
#include <QStringList>
#include <QVector>

class UTILSLIB_EXPORTS StringUtils  : public QObject {

    Q_OBJECT

public:

    static void removeUnwantedChars(const QString &allowedChars, QString *stringToClean);


    static QVector<int> findIndexesOfCharacterInString(const QString &str, const QChar &chr);


    template<typename QEnumClass>
    static QString enumClassToString (const QEnumClass value) {
        const int valueInt = static_cast<typename std::underlying_type<QEnumClass>::type>(value);
        return QString(QMetaEnum::fromType<QEnumClass>().valueToKey(valueInt));
    }


    template<typename QEnum>
    static QString enumToString (const QEnum value) {
        return QString(QMetaEnum::fromType<QEnum>().valueToKey(value));
    }


    static bool stringsMatch(const QString &s1, const QString &s2, bool caseSensitive);


};

#endif // STRINGUTILS_H
