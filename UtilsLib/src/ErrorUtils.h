// Created by anichols on 11/07/2021.

#ifndef ERRORUTILS_H
#define ERRORUTILS_H

#include "UtilsLib_Exports.h"

#include "Error.h"
#include "MathUtils.h"

using namespace Error;

enum class UTILSLIB_EXPORTS ErrorUtilsParam {
    ExcludeThreshold
    , IncludeThreshold
};


//TODO update ever method here to include optional error message.
class UTILSLIB_EXPORTS ErrorUtils {

public:

    /*!
    * @brief attempts to convert object to double.  Will return specified error if not possible
    */
    template<typename T>
    static Err toDouble(const T &object, double *out, Err e = eError) {
        bool ok;
        *out = object.toDouble(&ok);
        if (!ok) {
            qDebug() << object;
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief attempts to convert object to integer.  Will return specified error if not possible
    */
    template<typename T>
    static Err toInt(const T &object, int *out, Err e = eError) {
        bool ok;
        *out = object.toInt(&ok);
        if (!ok) {
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief Checks if a container is empty. Will return specified error if container is empty.
    */
    template <typename T>
    static Err isNotEmpty(const T &container, Err e = eEmptyContainerError) {
        if (container.isEmpty()) {
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief returns specified error if values are not equal
    */
    template<typename T>
    static Err isEqual(T s1, T s2, Err e = eError) {
        if (!MathUtils::tSame(s1, s2)) {
            qCritical()  << "Values are not equal" << s1 << "vs" << s2;
            rrr(e);
        }
        return eNoError;
    }

    template<typename T>
    static Err isEqualString(const T &s1, const T &s2, Err e = eError) {
        if(s1 != s2){
            return e;
        }

        return eNoError;
    }

    /*!
    * @brief returns specified error if value is false
    */
    static Err isTrue(bool s1, Err e = eError) {
        return isEqual(s1, true, e);
    }

    /*!
    * @brief returns specified error if value is false
    */
    static Err isFalse(bool s1, Err e = eError) {
        return isEqual(s1, false, e);
    }

    template<typename T>
    static Err isNotEqualString(const T &s1, const T &s2, Err e = eError) {
        if(s1 == s2){
            return e;
        }

        return eNoError;
    }

    /*!
    * @brief returns specified if values are equal
    */
    template<typename T>
    static Err isNotEqual(T s1, T s2, Err e = eError) {
        if (s1 == s2) {
            qCritical() << s1 << s2;
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief returns specified if index is bigger than size of container
    */
    template <typename T>
    static Err isIndexInRange(const T &checkThisVector, int index, Err e = eError) {
        if (index < 0 || index >= checkThisVector.size()) {
            qCritical() << index << checkThisVector.size();
            rrr(e);
        }

        return eNoError;
    }

    /*!
     * @brief returns specified if a given value is not between a specified start stop range.
     * Set includeBounds to true to include start and stop.
     */
    template <typename T>
    static Err isWithinRange(const T &value,
                             const T &start,
                             const T &end,
                             ErrorUtilsParam includeThreshold,
                             Err e = eError) {

        if (includeThreshold == ErrorUtilsParam::IncludeThreshold) {
            if (!(start < value && value < end) && !MathUtils::tSame(value, start)
            && !MathUtils::tSame(value, end)) {
                qCritical() << start << value << end;
                rrr(e);
            }
        } else {
            if (!(start < value && value < end)) {
                qCritical() << start << value << end;
                rrr(e);
            }
        }

        return eNoError;
    }

    /*!
    * @brief returns specified error if a given value is below threshold
    * Set includeThreshold to true to include threshold value in check
    */
    template <typename T>
    static Err isAboveThreshold(const T &value,
                                const T &threshold,
                                ErrorUtilsParam includeThreshold,
                                Err e = eError) {

        if (includeThreshold == ErrorUtilsParam::IncludeThreshold) {
            if (value < threshold) {
                qCritical() << value << threshold;
                rrr(e);
            }
        } else {
            if (value < threshold || MathUtils::tSame(value, threshold)) {
                qCritical() << value << threshold;
                rrr(e);
            }
        }

        return eNoError;
    }


    /*!
    * @brief returns specified error if a given value is above threshold
    * Set includeThreshold to true to include threshold value in check
    */
    template <typename T>
    static Err isBelowThreshold(const T &value,
                                const T &threshold,
                                ErrorUtilsParam includeThreshold,
                                Err e = eError) {

        if (includeThreshold == ErrorUtilsParam::IncludeThreshold) {
            if (value > threshold) {
                qCritical() << "Error vals for isBelowThreshold" << value << threshold;
                rrr(e);
            }
        } else {
            if (value > threshold || MathUtils::tSame(value, threshold)) {
                qCritical() << "Error vals for isBelowThreshold" << value << threshold;
                rrr(e);
            }
        }

        return eNoError;
    }


    static Err fileExists(const QString &filePath);


    template<typename T> //TODO make unit test
    static Err contains(T value, const QVector<T> &container, Err e = eError) {

        if (!container.contains(value)) {
            qCritical() << "Value is not in container";
            rrr(e);
        }

        return eNoError;
    }

};

#endif // ERRORUTILS_H
