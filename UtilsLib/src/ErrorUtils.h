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


class UTILSLIB_EXPORTS ErrorUtils {

public:

    /*!
    * @brief  Converts an object of type T to a double.
    * @tparam T: The datatype of the object. The object should support toDouble() conversion method.
    * @param object: The object that is to be converted to double.
    * @param out: Pointer to a double where the result of conversion will be stored.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after attempt of operation. Error could occur if conversion is not possible.
    * This error status should follow your internally defined error schemas/codes
    */
    template<typename T>
    static Err toDouble(const T &object, double *out, Err e = eError) {
        bool ok;
        *out = object.toDouble(&ok);
        if (!ok) {
            qDebug() << "input to convert" << object;
            rrr(e);
        }
        return eNoError;
    }

    template<typename T>
    static Err toLong(const T &object, long *out, Err e = eError) {
        bool ok;
        *out = object.toLong(&ok);
        if (!ok) {
            qDebug() << "input to convert" << object;
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Converts an object of type T to a float.
    * @tparam T: The datatype of the object. The object should support toDouble() conversion method.
    * @param object: The object that is to be converted to double.
    * @param out: Pointer to a double where the result of conversion will be stored.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after attempt of operation. Error could occur if conversion is not possible.
    * This error status should follow your internally defined error schemas/codes
    */
    template<typename T>
    static Err toFloat(const T &object, float *out, Err e = eError) {
        bool ok;
        *out = object.toFloat(&ok);
        if (!ok) {
            qDebug() << "input to convert" << object;
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Converts an object of type T to an integer.
    * @tparam T: The datatype of the object. The object should support toInt() conversion method.
    * @param object: The object that is to be converted to an integer.
    * @param out: Pointer to an int where the result of conversion will be stored.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after attempt of operation. Error could occur if conversion is not possible.
    * This error status should follow your internally defined error schemas/codes
    */
    template<typename T>
    static Err toInt(const T &object, int *out, Err e = eError) {
        bool ok;
        *out = object.toInt(&ok);
        if (!ok) {
            qDebug() << "input to convert" << object;
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Checks if a QVector of type T is not empty.
    * @tparam T: The datatype of the elements in QVector.
    * @param container: The QVector that is to be checked for emptiness.
    * @param e: Initial error status, eEmptyContainerError by default.
    * @return Err: Error status after the check. Error status will be 'eEmptyContainerError' if the container is empty and 'eNoError' if not empty.
    */
    template <typename T>
    static Err isNotEmpty(const QVector<T> &container, Err e = eEmptyContainerError) {
        if (container.isEmpty()) {
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type T is not empty.
    * @tparam T: The type of the container. The container should have an isEmpty() method for the check.
    * @param container: The container that is to be checked for emptiness.
    * @param e: Initial error status, eEmptyContainerError by default.
    * @return Err: Error status after the check. Error status will be 'eEmptyContainerError' if the container is empty and 'eNoError' if not empty.
    */
    template <typename T>
    static Err isNotEmpty(const T &container, Err e = eEmptyContainerError) {
        if (container.isEmpty()) {
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Checks if a std::vector of type T is not empty.
    * @tparam T: The datatype of the elements in std::vector.
    * @param container: The std::vector that is to be checked for emptiness.
    * @param e: Initial error status, eEmptyContainerError by default.
    * @return Err: Error status after the check. Error status will be 'eEmptyContainerError' if the container is empty and 'eNoError' if not empty.
    */
    template <typename T>
    static Err isNotEmpty(const std::vector<T> &container, Err e = eEmptyContainerError) {
        if (container.empty()) {
            rrr(e);
        }
        return eNoError;
    }

    /*!
    * @brief  Checks if two values of type T are equal.
    * @tparam T: The datatype of the values. The datatype T should be comparable using the MathUtils::tSame() function.
    * @param s1: The first value.
    * @param s2: The second value.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the values are not equal and 'eNoError' if they are equal.
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
	static Err isByteAligned(
		T* arr,
		size_t byteSize
		) {

    	ERR_INIT

	   e = isFalse(
		   reinterpret_cast<uintptr_t>(arr) % byteSize != 0,
		   eAVXAlignmentError
		   ); ree;

    	ERR_RETURN
	}

    /*!
    * @brief  Checks if two string-type values are equal.
    * @tparam T: The type of the string values. The type T should be comparable using the '!=' operator.
    * @param s1: The first string value.
    * @param s2: The second string value.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the strings are not equal and 'eNoError' if they are equal.
    */
    template<typename T>
    static Err isEqualString(const T &s1, const T &s2, Err e = eError) {
        if(s1 != s2){
            return e;
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a boolean value is true.
    * @param s1: The boolean value to check.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the value is false and 'eNoError' if it is true.
    * The check is done using the isEqual() function.
    */
    static Err isTrue(bool s1, Err e = eError) {
        return isEqual(s1, true, e);
    }

    /*!
    * @brief  Checks if a boolean value is false.
    * @param s1: The boolean value to check.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the value is true and 'eNoError' if it is false.
    * The check is done using the isEqual() function.
    */
    static Err isFalse(bool s1, Err e = eError) {
        return isEqual(s1, false, e);
    }

    /*!
    * @brief  Checks if two values of type T are not equal.
    * @tparam T: The datatype of the values. The datatype T should be comparable using the '==' operator.
    * @param s1: The first value.
    * @param s2: The second value.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the values are equal and 'eNoError' if they are not equal.
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
    * @brief  Checks if an index is in range for a given container of type T.
    * @tparam T: The type of the container. The container should have a size() method for the check.
    * @param checkThisVector: The container for which the index range has to be checked.
    * @param index: The index to be checked.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the index is not in range and 'eNoError' if it is in range.
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
    * @brief  Checks if a given value is within a specified range.
    * @tparam T: The datatype of the value, start, and end. The datatype should be suitable for comparison operations and should support the MathUtils::tSame() function if 'includeThreshold' is 'IncludeThreshold'.
    * @param value: The value to check whether it is within range.
    * @param start: The start of the range.
    * @param end: The end of the range.
    * @param includeThreshold: Specifies whether the 'start' and 'end' values are included in the range.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the value is not within the specified range and 'eNoError' if it is within the range.
    */
    template <typename T>
    static Err isWithinRange(
            const T &value,
            const T &start,
            const T &end,
            ErrorUtilsParam includeThreshold,
            Err e = eError
                    ) {

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

    template <typename T, typename U>
    static Err isWithinRange(
            const T &value,
            std::vector<U> &vec,
            Err e = eError
            ){

            if (value < 0 || value >= vec.size()) {
                qCritical() << value << "vs." << vec.size();
                rrr(e);
            }

            return eNoError;
        }

    /*!
    * @brief  Checks if a given value is above a specified threshold.
    * @tparam T: The datatype of the value and the threshold. The datatype should be suitable for comparison operations and should support the MathUtils::tSame() function if 'includeThreshold' is set to 'IncludeThreshold'.
    * @param value: The value to check whether it is above the threshold.
    * @param threshold: The specified threshold.
    * @param includeThreshold: Specifies whether the threshold is considered as 'above' (i.e., if the value is equal to the threshold, it's considered as being above the threshold).
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the value is not above the specified threshold and 'eNoError' if it is above the threshold.
    */
    template <typename T>
    static Err isAboveThreshold(
            const T &value,
            const T &threshold,
            ErrorUtilsParam includeThreshold,
            Err e = eError
                    ) {

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
    * @brief  Checks if a given value is below a specified threshold.
    * @tparam T: The datatype of the value and the threshold. The datatype should be suitable for comparison operations and should support the MathUtils::tSame() function if 'includeThreshold' is set to 'IncludeThreshold'.
    * @param value: The value to check whether it is below the threshold.
    * @param threshold: The specified threshold.
    * @param includeThreshold: Specifies whether the threshold is considered as 'below' (i.e., if the value is equal to the threshold, it's considered as being below the threshold).
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the value is not below the specified threshold and 'eNoError' if it is below the threshold.
    */
    template <typename T>
    static Err isBelowThreshold(
            const T &value,
            const T &threshold,
            ErrorUtilsParam includeThreshold,
            Err e = eError
                    ) {

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

    /*!
    * @brief Checks if a file exists at a given file path.
    * @param filePath: The file path to check.
    * @return Err: Error status after the check. Error status will be 'eFileError' if the file does not exist and 'eNoError' if it does exist.
    * The function first checks if the filePath is a non-empty string, and then checks if a file exists at the given filePath.
    */
    static Err fileExists(const QString &filePath);

    /*!
    * @brief  Checks if a container of type QVector<T> contains a specified value.
    * @tparam T: The datatype of the value and the elements of the QVector. The datatype should be suitable for comparison operations.
    * @param value: The value to check for its presence in the container.
    * @param container: The QVector container in which to search for the value.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the container does not contain the value and 'eNoError' if it does contain the value.
    */
    template<typename T>
    static Err contains(T value, const QVector<T> &container, Err e = eError) {

        if (!container.contains(value)) {
            qCritical() << "Value is not in container";
            rrr(e);
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type QMap<T, U> contains a specified key.
    * @tparam T: The datatype of the key and the keys in the QMap.
    * @tparam U: The datatype of the values in the QMap. This is not used in the function, but is required for specifying the QMap.
    * @param value: The key to check for its presence in the container.
    * @param container: The QMap container in which to search for the key.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the QMap does not contain the key and 'eNoError' if it does contain the key.
    */
    template<typename T, typename U>
    static Err contains(T value, const QMap<T, U> &container, Err e = eError) {

        if (!container.contains(value)) {
            qCritical() << "Value" << value << "is not in container";
            rrr(e);
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type QHash<T, U> contains a specified key.
    * @tparam T: The datatype of the key and the keys in the QHash.
    * @tparam U: The datatype of the values in the QHash. This is not used in the function, but is required for specifying the QHash.
    * @param value: The key to check for its presence in the container.
    * @param container: The QHash container in which to search for the key.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the QHash does not contain the key and 'eNoError' if it does contain the key.
    */
    template<typename T, typename U>
    static Err contains(T value, const QHash<T, U> &container, Err e = eError) {

        if (!container.contains(value)) {
            qCritical() << "Value" << value << "is not in container";
            rrr(e);
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type QVector<T> does not contain a specified value.
    * @tparam T: The datatype of the value and the elements in the QVector. The datatype should be suitable for comparison operations.
    * @param value: The value to check for its absence in the container.
    * @param container: The QVector container in which to check for the absence of the value.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the container does contain the value and 'eNoError' if it does not contain the value.
    */
    template<typename T>
    static Err doesNotContain(T value, const QVector<T> &container, Err e = eError) {

        if (container.contains(value)) {
            qCritical() << "Value is not in container";
            rrr(e);
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type QMap<T, U> does not contain a specified key.
    * @tparam T: The datatype of the key and the keys in the QMap.
    * @tparam U: The datatype of the values in the QMap. This is not used in the function, but is required for specifying the QMap.
    * @param value: The key to check for its absence in the container.
    * @param container: The QMap container in which to check for the absence of the key.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the QMap does contain the key and 'eNoError' if it does not contain the key.
    */
    template<typename T, typename U>
    static Err doesNotContain(T value, const QMap<T, U> &container, Err e = eError) {

        if (container.contains(value)) {
            qCritical() << "Value is not in container";
            rrr(e);
        }

        return eNoError;
    }

    /*!
    * @brief  Checks if a container of type QHash<T, U> does not contain a specified key.
    * @tparam T: The datatype of the key and the keys in the QHash.
    * @tparam U: The datatype of the values in the QHash. This is not used in the function but is required for specifying the QHash.
    * @param value: The key to check for its absence in the container.
    * @param container: The QHash container in which to check for the absence of the key.
    * @param e: Initial error status, eError by default.
    * @return Err: Error status after the check. Error status will be 'eError' if the QHash does contain the key and 'eNoError' if it does not contain the key.
    */
    template<typename T, typename U>
    static Err doesNotContain(T value, const QHash<T, U> &container, Err e = eError) {

        if (container.contains(value)) {
            qCritical() << "Value is not in container";
            rrr(e);
        }

        return eNoError;
    }

};

#endif // ERRORUTILS_H
