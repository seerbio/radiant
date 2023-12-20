//
// Created by anichols on 12/8/23.
//

#ifndef PYTHIADIACPP_POINTFF_H
#define PYTHIADIACPP_POINTFF_H

#include "UtilsLib_Exports.h"
#include "MathUtils.h"

#include <QDataStream>
#include <QDebug>

class UTILSLIB_EXPORTS PointFF {

public:

    PointFF();
    PointFF(float x, float y);

    ~PointFF() = default;

    float x();
    float y();

    [[nodiscard]] float x() const;
    [[nodiscard]] float y() const;

    float& rx();
    float& ry();

    friend QDebug operator<<(QDebug dbg, const PointFF& obj) {
        dbg.nospace() << "PointFF(" << obj.x() << ", " << obj.y() << ") ";
        return dbg;
    }

    friend QDataStream& operator<<(QDataStream &stream, const PointFF &obj) {
        stream << obj.x() << obj.y();
        return stream;
    }

    friend QDataStream& operator>>(QDataStream &stream, PointFF &obj) {
        stream >> obj.rx() >> obj.ry();
        return stream;
    }

    bool operator==(const PointFF& other) const {

        if (MathUtils::tSame(x(), other.x()) && MathUtils::tSame(y(), other.y())) {
            return true;
        }

        return false;
    }

    bool operator!=(const PointFF& other) const {

        if (*this == other) {
            return false;
        }

        return true;
    }

private:

    float m_x;
    float m_y;

};


#endif //PYTHIADIACPP_POINTDF_H
