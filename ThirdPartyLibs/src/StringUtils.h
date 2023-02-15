#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <iostream>
#include <vector>


class  StringUtils {

public:

    static std::vector<std::string> split(
            const std::string &str,
            const std::string &delim
            );

    static std::string replaceSubstring(
            std::string str,
            const std::string &to_replace,
            const std::string &replace_with
            );

//    static void removeUnwantedChars(
//            const std::string &allowedChars,
//            std::string *stringToClean
//    );
//
//
//    static std::vector<int> findIndexesOfCharacterInString(
//            const QString &str,
//            const QChar &chr
//    );





//    static bool stringsMatch(
//            const std::string &s1,
//            const std::string &s2,
//            bool caseSensitive
//    );

};

#endif // STRINGUTILS_H
