#include "StringUtils.h"


void StringUtils::removeUnwantedChars(const QString &allowedChars, QString *stringToClean)
{

    QStringList stringToCleanList = stringToClean->split("");

    auto deteteIters = std::remove_if(stringToCleanList.begin(), stringToCleanList.end(),
                                      [allowedChars](const QString &c){return !allowedChars.contains(c);});

    stringToCleanList.erase(deteteIters, stringToCleanList.end());
    *stringToClean = stringToCleanList.join("");
}


QVector<int> StringUtils::findIndexesOfCharacterInString(const QString &str, const QChar &chr) {

    QVector<int> indexes;

    for (int i = 0; i < str.size(); i++) {

        const QChar &c = str.at(i);
        if (c != chr) {
            continue;
        }

        indexes.push_back(i);
    }

    return indexes;
}


bool StringUtils::stringsMatch(const QString &s1, const QString &s2, bool caseSensitive) {
    Qt::CaseSensitivity sensitivity = caseSensitive ? Qt::CaseSensitive : Qt::CaseInsensitive;
    return QString::compare(s1, s2,  sensitivity) == 0;
}
