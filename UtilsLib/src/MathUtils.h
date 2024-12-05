//
// Created by anichols on 11/07/2021.
//

#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "UtilsLib_Exports.h"

#include "Error.h"

#include <QDebug>
#include <QElapsedTimer>
#include <QMap>
#include <QPointF>

#include <cmath>
#include <vector>

//NOTE: ErrorUtils is not used here because something funky happens when it's imported.

using namespace Error;

namespace MathUtilsConstants {

    extern const double UTILSLIB_EXPORTS PI;


}//NAMESPACE


class UTILSLIB_EXPORTS MathUtils {


public:

    /*!
    * @brief  Calculates the median value of a given numeric container.
    * @tparam Vector: The type of the numeric container. The container should support size(), indexing, and sorting operations.
    * @param vecOriginal: The numeric container for which to calculate the median.
    * @return double: The median of the given container.
    * If the number of elements in the vector is 0, the function will return default double value (0.0).
    * If the number of elements in the vector is 1 or 2, the function will return the reasonable statistic.
    * If the number of elements in the vector is 3, the function will return the value of the middle element.
    * If the number of elements in the vector is more than 3, the function will either return the value of the middle element (if the count of elements in vector is even)
    * or average of two middle elements (if the count of elements in vector is odd).
    */
    template <typename Vector>
    static double median(const Vector &vecOriginal) {
        Vector vec = vecOriginal;
        const auto vecSize = static_cast<size_t>(vec.size());

        switch (vecSize) {
            case 0:
                return double();
            case 1:
                return vec[0];
            case 2:
                return (vec[0] + vec[1]) * 0.5;
            case 3:
                std::sort(vec.begin(), vec.end());
                return vec[1];
            default:;
        }

        const int buffer = 1;
        const size_t midPoint = vecSize / 2;

        std::partial_sort(vec.begin(),
                          vec.begin() + midPoint + buffer,
                          vec.end());

        if (vecSize % 2 != 0) {
            return vec[midPoint];
        }

        return (vec[midPoint] + vec[midPoint - 1]) * 0.5;
    }

    /*!
    * @brief  Calculates the mean (average) of a given numeric container.
    * @tparam Vector: The type of numeric container. The container should support size(), begin(), end() operations
    *      and the elements should be numeric (supporting addition and division).
    * @param vec: The numeric container for which to calculate the mean.
    * @return double: The mean of the given container. If the vector is empty, the function will return 0.0.
    */
    template <typename Vector>
    static double mean(const Vector &vec) {
        if(vec.empty()){
            return 0.0;
        }

        return std::accumulate(vec.begin(), vec.end(), 0.0) / static_cast<double>(vec.size());
    }

    /*!
     * @brief  Calculates the weighted mean of a provided list of values and weights.
     * @tparam T: The datatype of the values and weights. The datatype should be numeric (int, float, double, etc.) and support arithmetic operations.
     * @param vecVals: A standard vector representing the numeric values.
     * @param vecWeights: A standard vector representing the weights of the numeric values.
     * @param inverseWeights: A boolean flag that determines whether to inverse the weights or not.
     * @param weightedAverage: A pointer to a double which the function will use to output the calculated weighted mean.
     * @return Err: Error status after the checks and calculations. Error status will be 'eValueError' if the value and weight vectors are of different sizes, sum of weights is zero, or 'eNoError' otherwise.
     */
    template <typename T>
    static Err weightedMean(
            const std::vector<T> &vecVals,
            const std::vector<T> &vecWeights,
            bool inverseWeights,
            double *weightedAverage
            ) {

        ERR_INIT

        if (vecVals.size() != vecWeights.size()) {
            rrr(eValueError);
        }

        double sum = 0;
        double weightSum = 0;

        for(int i = 0; i < vecVals.size(); i++) {

            T w = std::max(std::abs(vecWeights.at(i)), 0.01); //TODO figure out a better solution.

            const double weightVal = inverseWeights ? (1.0 / w) : w;
            sum += vecVals.at(i) * weightVal;
            weightSum += vecWeights.at(i);
        }

        if (MathUtils::tZero(weightSum)) {
            rrr(eValueError);
        }

        *weightedAverage = sum / weightSum;
        ERR_RETURN
    }

    /*!
    * @brief  Calculates the standard deviation of a given numeric container.
    * @tparam Vector: The type of the numeric container. The container should support size(), begin(), end() operations and the elements should be numeric (supporting addition, subtraction, squaring, and square root operations).
    * @param vec: The numeric container for which to calculate the standard deviation.
    * @return double: The standard deviation of the given container. The function first calculates the mean of the elements, then it computes the square of the difference between each element and the mean.
    * Finally, it calculates the square root of the average of these squared differences.
    */
    template <typename Vector>
    static double stDev(const Vector &vec) {

        const double vecMean = mean(vec);

        const auto transformLogic = [vecMean](double vecVal) {
            const double difference = vecVal - vecMean;
            return  pow(difference, 2);
        };

        QVector<double> differences;
        std::transform(vec.begin(), vec.end(), std::back_inserter(differences), transformLogic);

        const double differencesSum = std::accumulate(differences.begin(),  differences.end(), 0.0);

        return std::sqrt(differencesSum / vec.size());
    }

    /*!
    * @brief  Checks if value is very close to zero (zero in terms of fuzzy comparison)
    * @tparam T: The datatype of the value that needs to be checked. The datatype should be numeric (float, double, etc.) that can be compared with zero.
    * @param val: The numeric value that needs to be checked.
    * @return bool: Returns true if the value is very close to zero; false otherwise.
    *
    * Note: The comparison is done using the qFuzzyIsNull() function from Qt which considers very small floating-point values near zero as zero.
    */
    template <typename T>
    static bool tZero(T val) {
        return qFuzzyIsNull(val);
    }

    /*!
    * @brief  Checks if two numbers are "approximately equal" or "same" within a small difference (fudgeFactor).
    * @tparam T: The datatype of the values to be checked. The datatype should be numeric that can support subtraction and comparison with another numeric type (fudgeFactor).
    * @param v1: The first value to be compared.
    * @param v2: The second value to be compared.
    * @param fudgeFactor: The tolerance value within which v1 and v2 are considered the same. Defaults to 0.01 if not specified.
    * @return bool: Returns true if the absolute difference between v1 and v2 is less than fudgeFactor; false otherwise.
    *
    * Note: The comparison takes into account a small tolerance and two numbers are considered "same" if they are close enough.
    */
    template <typename T>
    static bool tSame(T v1, T v2, double fudgeFactor = 0.01) {
        //TODO Figure out a better way to do this.
        return std::abs(v1 - v2) < fudgeFactor;
    }

    template <typename T>
    static T pRound(T val, int precision = 1) {
        const T precisionMultiplier = pow(10, precision);
        return round(val * precisionMultiplier) / precisionMultiplier;
    }

    /*!
    * @brief  Rounds a number to a specified number of decimal places.
    * @tparam T: The datatype of the value to be rounded. The datatype should be numeric (int, float, double, etc.) and should support arithmetic operations.
    * @param val: The numeric value that needs to be rounded.
    * @param precision: The number of decimal places to which to round. Defaults to 1 if not specified.
    * @return T: The rounded numeric value.
    *
    * Note: The rounding is done by multiplying the value with the precision factor (10^precision), then rounding to the nearest integer, and finally dividing by the precision factor.
    */
    template<typename T>
    static int hashDecimal(T val, int precision = 1) {
        const double precisionMultiplier = pow(10, precision);
        return static_cast<int>(round(val * precisionMultiplier));
    }

    /*!
    * @brief  Hashes a decimal number to an integer using a specified precision.
    * @tparam T: The datatype of the value to be hashed. The datatype should be numeric (float, double, etc.) and support arithmetic operations.
    * @param val: The decimal value that needs to be hashed.
    * @param precision: The precision to be used for hashing. Provides the level of granularity of the hash.
    * @return int: The hashed value as an integer.
    *
    * Note: The value is hashed by multiplying it with the reciprocal of the precision, rounding it to the nearest integer, and then casting it to an integer.
    */
    template<typename T>
    static int hashDecimal(T val, double precision) {
        const int precisionMultiplier = static_cast<int>(round(1.0 / precision));
        return static_cast<int>(round(val * precisionMultiplier));
    }

    template<typename T>
    static T unHashDecimal(int val, int precision = 1) {
        const double precisionMultiplier = pow(10, precision);
        return static_cast<T>(val / precisionMultiplier);
    }

    /*!
    * @brief  Converts a hashed decimal number back to its original form using a specified precision.
    * @tparam T: The datatype of the original decimal value before hashing. The datatype should be numeric (float, double, etc.) and support arithmetic operations.
    * @param val: The hashed integer value that needs to be converted back.
    * @param precision: The precision that was used during hashing. Defaults to 1 if not specified.
    * @return T: The original decimal number converted back from the hash.
    *
    * Note: The conversion is done by dividing the hashed value by the precision factor (10^precision) and then casting it back to the original datatype.
    */
    template<typename T>
    static T calculatePPM(T val, T ppmTolerance) {
        return (val * ppmTolerance) / 1e6;
    }

    /*!
    * @brief Calculates the factorial of a given number.
    * @param n: The positive integer whose factorial is to be calculated.
    * @return double: The factorial of the input number. If the number is negative, the function will return 0.
    * If the calculated factorial is less than 0 (which can occur due to integer overflow for large input numbers),
    * the function will return the maximum value that a double can represent.
    * The factorial calculation is done using the recursive definition of factorial, i.e., n! = n*(n-1)!.
    */
    static double factorial(int n);

    /*!
    * @brief  Finds the index of the element in the container that is closest to a given value.
    * @tparam T: The datatype of the elements in the container. The datatype should be numeric (int, float, double, etc.) and it should support subtraction and absolute value operations.
    * @param vec: The numeric container (QVector) in which to find the closest element.
    * @param value: The value to which the closest element should be found.
    * @return int: Returns the index of the element that is closest to value. If the vector is empty, the function will return -1.
    *
    */
    template<typename T>
    static int closest(const QVector<T> &vec, T value) {

        assert(!vec.isEmpty());
        auto it = std::min_element(vec.begin(), vec.end(), [value] (T a, T b) {
            return std::abs(a - value) <= std::abs(b - value);
        });

        if(it == vec.end()) {
            return -1;
        }

        return it - vec.begin();
    }

    /*!
    * @brief  Finds the index of the element in the container that is closest to a given value.
    * @tparam T: The datatype of the elements in the container. The datatype should be numeric (int, float, double, etc.) and it should support subtraction and absolute value operations.
    * @param vec: The numeric container (std::vector) in which to find the closest element.
    * @param value: The value to which the closest element should be found.
    * @return int: Returns the index of the element that is closest to value. If the vector is empty, the function will return -1.
    *
    */
    template<typename T>
    static int closest(const std::vector<T> &vec, T value) {

        assert(!vec.isEmpty());
        auto it = std::min_element(vec.begin(), vec.end(), [value] (T a, T b) {
            return std::abs(a - value) <= std::abs(b - value);
        });

        if(it == vec.end()) {
            return -1;
        }

        return it - vec.begin();
    }

    /*!
    * @brief  Finds the QPointF in the container that has the "x" value closest to a given value.
    * @param vec: The QVector of QPoints in which to find the closest "x" value.
    * @param value: The "x" value to which the closest QPointF should be found.
    * @return QPointF: Returns the QPointF that has the "x" value closest to the input value.
    * If the vector is empty, the function will return a QPoint {-1, -1}.
    *
    */
    static QPointF closestXValPoint(const QVector<QPointF> &vec, double value);

    /*!
    * @brief  Finds the index of the maximum value in a given numeric container.
    * @tparam Vector: The type of the numeric container. The container should support size(), begin(), end() operations and the elements should be numeric (supporting addition and comparision operations).
    * @param vec: The numeric container for which to find the index of the maximum value.
    * @return int: The index of the maximum value in the container. If the sum of all elements in vector is approximately zero (checked by tZero function), or vector is empty, the function will return -1.
    *
    * Note: tZero function is used to compare the sum with zero which can account for very small floating point values and reduce them to zero.
    */
    template<typename Vector>
    static int findMaxIndexInVector(const Vector &vec) {

        const double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
        if (MathUtils::tZero(sum)) {
            return -1;
        }

        return std::max_element(vec.begin(), vec.end()) - vec.begin();
    }

    /*!
    * @brief  Calculates the root mean square error (RMSE) for a set of actual and predicted values.
    * @tparam T: The datatype of the actual and predicted values. The datatype should be numeric (supporting subtraction, squaring, division, and square root operations).
    * @param actualVsPredicted: A container of QPair elements where each pair contains actual and predicted values.
    * @return T: The root mean square error of the given values.
    * The function calculates the square of the difference between actual and predicted values, calculates the mean of these squared differences,
    * and finally returns the square root of this mean value.
    */
    template<typename T>
    static T rmse(const QVector<QPair<T, T>> &actualVsPredicted) {

        const auto accLogic = [](T sum, const QPair<T, T> &pr){
            return sum + std::pow((pr.first - pr.second), 2);
        };

        const double squaredDiffs = std::accumulate(
                actualVsPredicted.begin(),
                actualVsPredicted.end(),
                0.0,
                accLogic
                );

        return std::sqrt(squaredDiffs / actualVsPredicted.size());
    }

    /*!
    * @brief  Calculates qValue and decoyRatio for given identifiers.
    * @tparam Identifier: The datatype of the keys in the QHash containers, typically a string or an int.
    * @tparam T: The datatype of the values in the QHash containers. It should support the '<' operator, arithmetic operations and conversion to 'double'.
    * @param identifierVsTarget: A QHash container mapping identifiers to target values.
    * @param identifierVsDecoys: A QHash container mapping identifiers to decoy values.
    * @param identifierVsQValue: A pointer to a QHash container where mappings from identifiers to calculated qValues will be stored.
    * @param identifierVsDecoyRatio: A pointer to a QHash container where mappings from identifiers to calculated decoyRatios will be stored.
    * @return Err: Error status after the checks and calculations. Error status will be 'eValueError' if any of the input containers are empty and 'eNoError' otherwise.
    */
    template<typename Identifier, typename T>
    static Err calculateQValue(
            const QVector<QPair<Identifier, T>> &identifierVsTarget,
            const QHash<Identifier, T> &identifierVsDecoys,
            QHash<Identifier, QPair<T, T>> *identifierVsQValueVsDecoyRatio
            ) {

        ERR_INIT

        if (identifierVsTarget.isEmpty() || identifierVsDecoys.isEmpty()) {
            rrr(eValueError);
        }

        identifierVsQValueVsDecoyRatio->clear();

        QElapsedTimer et;
        et.start();

        QVector<T> targetScores;
        std::transform(
            identifierVsTarget.begin(),
            identifierVsTarget.end(),
            std::back_inserter(targetScores),
            [](const QPair<Identifier, T> &pr){return pr.second;}
            );
        std::sort(targetScores.begin(), targetScores.end());

        QVector<T> decoyScores = identifierVsDecoys.values().toVector();
        std::sort(decoyScores.begin(), decoyScores.end());

        for (const QPair<Identifier, T> &pr : identifierVsTarget) {

            const Identifier &identifier = pr.first;
            T score = pr.second;

            auto decoyIndexLowest = std::lower_bound(decoyScores.begin(), decoyScores.end(), score);
            if (decoyIndexLowest > decoyScores.begin()) {
                if (*(decoyIndexLowest - 1) > decoyScores.front()) {
                    score = *(--decoyIndexLowest);
                }
            }

            auto targetIndexLowest = std::lower_bound(targetScores.begin(), targetScores.end(), score);

            const long targetCount = std::distance(targetIndexLowest, targetScores.end());
            const long decoyCount = std::distance(decoyIndexLowest, decoyScores.end());

            const T qvalue = std::min(
                    1.0,
                    (static_cast<double>(std::max(static_cast<long>(0), decoyCount))) / static_cast<double>(std::max(static_cast<long>(1), targetCount))
                    );

            const T decoyRatio
                = std::min(1.0, (static_cast<double>(std::max(static_cast<long>(0), decoyCount)) / static_cast<double>(std::max(1, identifierVsTarget.size()))));
            identifierVsQValueVsDecoyRatio->insert(identifier, {qvalue, decoyRatio});

        }

        ERR_RETURN
    }

    /*!
    * @brief  Generates a list of unique random numbers and returns a QMap indicating which indices are included in the random selection.
    * @param totalSizeOfList: The total size of the list (upper bound of the range of random values).
    * @param desiredSizeRandomNumbers: The number of unique random numbers desired.
    * @param seed: The seed for the random number generator.
    * @return QMap<int, bool>: Returns a QMap where keys are positions from 0 to totalSizeOfList and values are boolean flags indicating whether the position is in the random selection (true) or not (false).
    *
    * Note: This method uses the Mersenne Twister engine (std::mt19937) to generate random numbers and std::uniform_int_distribution to distribute these numbers evenly within the specified range.
    */
    static QMap<int, bool> generateRandomSelectionList(
            int totalSizeOfList,
            int desiredSizeRandomNumbers,
            int seed = 666
            );

    /*!
    * @brief  Splits given data into training and test datasets according to a specified fraction.
    * @tparam T: The type of the data elements in the container. It can be any complex type that can be copied and pushed into QVector container.
    * @param data: The main dataset to be split, represented as a QVector of data elements.
    * @param testFraction: The fraction of the original dataset which should serve as the test dataset.
    * @param seed: The seed for the random generator when selecting test data.
    * @param trainingData: A pointer to a QVector where the generated training dataset will be stored.
    * @param testData: A pointer to a QVector where the generated test dataset will be stored.
    * @return Err: Error status after the checks and operations. Error status will be 'eValueError' if the data vector is empty or testFraction is zero, and 'eNoError' otherwise.
    *
    * Note: Training and test datasets are disjoint. Each data item will be included in either the training or the test dataset, but not both.
    */
    template<typename T>
    static Err trainTestSplit(
            const QVector<T> &data,
            double testFraction,
            int seed,
            QVector<T> *trainingData,
            QVector<T> *testData
            ) {

        ERR_INIT

        if (data.isEmpty() || MathUtils::tZero(testFraction)) {
            rrr(eValueError);
        }

        trainingData->clear();
        testData->clear();

        const int testDataSize = static_cast<int>(data.size() * testFraction);

        const QMap<int, bool> selectionList = MathUtils::generateRandomSelectionList(
                data.size(),
                testDataSize,
                seed
                );

        for (int i = 0; i < data.size(); i++) {

            const T &d = data.at(i);

            if (selectionList.value(i)) {
                testData->push_back(d);
                continue;
            }

            trainingData->push_back(d);
        }

        ERR_RETURN
    }

    /**
    * @brief Calculate mass accuracy parts per million (PPM)
    *
    * @tparam T Template type
    * @param massTheo Theoretical mass
    * @param massObserved Observed mass
    * @return Mass accuracy PPM calculated as 1e6 * (massObserved - massTheo) / massTheo
    */
    template<typename T>
    static T calculateMassAccuracyPPM(T massTheo, T massObserved) {
        return 1e6 * (massObserved - massTheo) / massTheo;
    }

    template<typename T>
    static bool vectorContainsInfOrNaN(const QVector<T>& vec) {
        for (T val : vec) {
            if (std::isinf(val) || std::isnan(val)) {
                return true;
            }
        }
        return false;
    }

};


#endif //MATHUTILS_H
