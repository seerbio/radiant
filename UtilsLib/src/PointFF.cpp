//
// Created by anichols on 12/8/23.
//

#include "PointFF.h"


PointFF::PointFF()
: m_x(0.0)
, m_y(0.0)
{}

PointFF::PointFF(float x, float y)
: m_x(x)
, m_y(y)
{}

float PointFF::x() {
    return m_x;
}

float PointFF::y() {
    return m_y;
}

float PointFF::x() const {
    return m_x;
}

float PointFF::y() const {
    return m_y;
}

float &PointFF::rx() {
    return m_x;
}

float &PointFF::ry() {
    return m_y;
}


