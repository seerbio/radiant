//
// Created by anichols on 12/8/23.
//

#include "PointDF.h"


PointDF::PointDF()
: m_x(0.0)
, m_y(0.0)
{}

PointDF::PointDF(double x, float y)
: m_x(x)
, m_y(y)
{}

double PointDF::x() {
    return m_x;
}

float PointDF::y() {
    return m_y;
}

double PointDF::x() const {
    return m_x;
}

float PointDF::y() const {
    return m_y;
}

double &PointDF::rx() {
    return m_x;
}

float &PointDF::ry() {
    return m_y;
}


