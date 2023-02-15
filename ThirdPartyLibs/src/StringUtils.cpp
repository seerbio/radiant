#include "StringUtils.h"

#include <sstream>

std::vector<std::string> StringUtils::split(const std::string &str, const std::string &delim) {

    std::vector<std::string> res;
    size_t pos = 0;

    while (pos < str.size()) {
        size_t found = str.find(delim, pos);
        if (found == std::string::npos) {
            res.push_back(str.substr(pos));
            break;
        }
        res.push_back(str.substr(pos, found - pos));
        pos = found + delim.size();
    }

    return res;
}


std::string StringUtils::replaceSubstring(
        std::string str,
        const std::string &to_replace,
        const std::string &replace_with
) {

    size_t start_pos = 0;

    while ((start_pos = str.find(to_replace, start_pos)) != std::string::npos) {
        str.replace(start_pos, to_replace.length(), replace_with);
        start_pos += replace_with.length();
    }

    return str;
}


//void StringUtils::removeUnwantedChars(
//        const std::string &allowedChars,
//        std::string *stringToClean
//) {
//
//    QStringList stringToCleanList = stringToClean->split("");
//
//    auto deteteIters = std::remove_if(stringToCleanList.begin(), stringToCleanList.end(),
//                                      [allowedChars](const QString &c){return !allowedChars.contains(c);});
//
//    stringToCleanList.erase(deteteIters, stringToCleanList.end());
//    *stringToClean = stringToCleanList.join("");
//}
//
//
//QVector<int> StringUtils::findIndexesOfCharacterInString(
//        const QString &str,
//        const QChar &chr
//) {
//
//    QVector<int> indexes;
//
//    for (int i = 0; i < str.size(); i++) {
//
//        const QChar &c = str.at(i);
//        if (c != chr) {
//            continue;
//        }
//
//        indexes.push_back(i);
//    }
//
//    return indexes;
//}
//
//
//bool StringUtils::stringsMatch(
//        const QString &s1,
//        const QString &s2,
//        bool caseSensitive
//) {
//
//    Qt::CaseSensitivity sensitivity = caseSensitive ? Qt::CaseSensitive : Qt::CaseInsensitive;
//
//    return QString::compare(s1, s2,  sensitivity) == 0;
//}
