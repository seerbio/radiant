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

using namespace Error;

namespace MathUtilsConstants {

    extern const double UTILSLIB_EXPORTS PI;


}//NAMESPACE


class UTILSLIB_EXPORTS MathUtils {


public:

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

    template <typename Vector>
    static double mean(const Vector &vec) {
        if(vec.empty()){
            return 0.0;
        }

        return std::accumulate(vec.begin(), vec.end(), 0.0) / static_cast<double>(vec.size());
    }

    //TODO write test for this
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

    template <typename T>
    static bool tZero(T val)
    {
        return qFuzzyIsNull(val);
    }

    template <typename T>
    static bool tSame(T v1, T v2, double fudgeFactor = 0.01)
    {
        //TODO Figure out a better way to do this.
        return std::abs(v1 - v2) < fudgeFactor;
    }

    template <typename T>
    static T pRound(T val, int precision = 1)
    {
        const double precisionMultiplier = pow(10, precision);
        return round(val * precisionMultiplier) / precisionMultiplier;
    }

    template<typename T>
    static int hashDecimal(T val, int precision = 1) {
        const double precisionMultiplier = pow(10, precision);
        return static_cast<int>(round(val * precisionMultiplier));
    }

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

    static double calculatePPM(double val, double ppmTolerance);

    static double factorial(int n);

    template<typename T>
    static int closest(const QVector<T> &vec, T value) {

        assert(!vec.isEmpty());
        auto it = std::min_element(vec.begin(), vec.end(), [value] (double a, double b) {
            return std::abs(a - value) <= std::abs(b - value);
        });

        if(it == vec.end()) {
            return -1;
        }

        return it - vec.begin();
    }

    static QPointF closestXValPoint(const QVector<QPointF> &vec, double value);

    template<typename T>
    static int findMaxIndexInVector(const QVector<T> &vec) {

        const int sum = std::accumulate(vec.begin(), vec.end(), 0);
        if (sum == 0) {
            return -1;
        }

        return std::max_element(vec.begin(), vec.end()) - vec.begin();
    }

    template<typename T>
    static double rmse(const QVector<QPair<T, T>> &actualVsPredicted) {

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

    template<typename Identifier>
    static Err calculateQValue(
            const QHash<Identifier, double> &identifierVsTarget,
            const QHash<Identifier, double> &identifierVsDecoys,
            QHash<Identifier, double> *identifierVsQValue,
            QHash<Identifier, double> *identifierVsDecoyRatio
            ) {

        ERR_INIT

        if (identifierVsTarget.isEmpty() || identifierVsDecoys.isEmpty()) {
            rrr(eValueError);
        }

        identifierVsQValue->clear();
        identifierVsDecoyRatio->clear();

        QElapsedTimer et;
        et.start();

        QVector<double> targetScores = identifierVsTarget.values().toVector();
        std::sort(targetScores.begin(), targetScores.end());

        QVector<double> decoyScores = identifierVsDecoys.values().toVector();
        std::sort(decoyScores.begin(), decoyScores.end());

        for (auto it = identifierVsTarget.begin(); it != identifierVsTarget.end(); it++) {

            const Identifier &identifier = it.key();
            double score = it.value();

            auto decoyIndexLowest = std::lower_bound(decoyScores.begin(), decoyScores.end(), score);
            if (decoyIndexLowest > decoyScores.begin()) {
                if (*(decoyIndexLowest - 1) > decoyScores.front()) {
                    score = *(--decoyIndexLowest);
                }
            }

            auto targetIndexLowest = std::lower_bound(targetScores.begin(), targetScores.end(), score);

            const long targetCount = std::distance(targetIndexLowest, targetScores.end());
            const long decoyCount = std::distance(decoyIndexLowest, decoyScores.end());

            const double qvalue
                = std::min(1.0, (static_cast<double>(std::max(static_cast<long>(1), decoyCount))) / static_cast<double>(std::max(static_cast<long>(1), targetCount)));
            identifierVsQValue->insert(identifier, qvalue);

            const double decoyRatio
                = std::min(1.0, (static_cast<double>(std::max(static_cast<long>(1), decoyCount)) / static_cast<double>(std::max(1, identifierVsTarget.size()))));
            identifierVsDecoyRatio->insert(identifier, decoyRatio);

        }

        ERR_RETURN
    }

    static QMap<int, bool> generateRandomSelectionList(
            int totalSizeOfList,
            int desiredSizeRandomNumbers,
            int seed = 666
            );

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

};


#endif //MATHUTILS_H
