//
// Created by anichols on 11/07/2021.
//

#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "UtilsLib_Exports.h"

#include <QDebug>

#include <cmath>
#include <vector>



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
    static double mean(const Vector &vec)
    {
        if(vec.empty()){
            return 0.0;
        }

        return std::accumulate(vec.begin(), vec.end(), 0.0) / static_cast<double>(vec.size());
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

    static unsigned long long factorial(int n);

    static double calculateHyperScore(long long intensity, int bIonCount, int yIonsCount);

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

};


#endif //MATHUTILS_H
