#ifndef TDOA_GRAD_H
#define TDOA_GRAD_H

#include "Point.h"

class grad {
public:
    // P = (x, y) аргумент
    // B C D -- заданные точки
    // AB_AC -- это AB - AC разность хода
    // AB_AD -- это AB - AD разность хода
    // F -- целевая функция
    static double F(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);

    // Частная производная dF/dx
    static double dF_dx(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);

    // Частная производная dF/dy
    static double dF_dy(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);

    // градиентный метод с дроблением шага
    static Point gradMethod(const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);
};


#endif //TDOA_GRAD_H
