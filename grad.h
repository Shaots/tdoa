#ifndef TDOA_GRAD_H
#define TDOA_GRAD_H

#include "point.h"

class grad {
public:
    // P = (x, y) аргумент
    // B C D -- заданные точки
    // AB_AC -- это AB - AC
    // AB_AD -- это AB - AD
    // F -- целевая функция
    static double F(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    // Частная производная dF/dx
    static double dF_dx(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    // Частная производная dF/dy
    static double dF_dy(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    // градиентный метод с дроблением шага
    static point gradMethod(const point &B, const point &C, const point &D, double AB_AC, double AB_AD);
};


#endif //TDOA_GRAD_H
