//
// Created by anichols on 12/8/23.
//

#ifndef PYTHIADIACPP_POINTDF_H
#define PYTHIADIACPP_POINTDF_H

#include "MathUtils.h"
#include "UtilsLib_Exports.h"

#include <QDataStream>
#include <QDebug>

class UTILSLIB_EXPORTS PointDF {

public:

    PointDF();
    PointDF(double x, float y);

    ~PointDF() = default;

    double x();
    float y();

    [[nodiscard]] double x() const;
    [[nodiscard]] float y() const;

    double& rx();
    float& ry();

    friend QDebug operator<<(QDebug dbg, const PointDF& obj) {
        dbg.nospace() << "PointDF(" << obj.x() << ", " << obj.y() << ") ";
        return dbg;
    }

    friend QDataStream& operator<<(QDataStream &stream, const PointDF &obj) {
        stream << obj.x() << obj.y();
        return stream;
    }

    friend QDataStream& operator>>(QDataStream &stream, PointDF &obj) {
        stream >> obj.rx() >> obj.ry();
        return stream;
    }

    bool operator==(const PointDF& other) const {

        if (MathUtils::tSame(x(), other.x()) && MathUtils::tSame(y(), other.y())) {
            return true;
        }

        return false;
    }

    bool operator!=(const PointDF& other) const {

        if (*this == other) {
            return false;
        }

        return true;
    }

private:

    double m_x;
    float m_y;

};


#endif //PYTHIADIACPP_POINTDF_H
