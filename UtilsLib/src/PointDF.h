//
// Created by anichols on 12/8/23.
//

#ifndef PYTHIADIACPP_POINTDF_H
#define PYTHIADIACPP_POINTDF_H


class PointDF {

public:

    PointDF();
    PointDF(double x, float y);

    ~PointDF() = default;

    double x();
    float y();

    double& rx();
    float& ry();

private:

    double m_x;
    float m_y;

};


#endif //PYTHIADIACPP_POINTDF_H
