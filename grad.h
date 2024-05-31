#ifndef TDOA_GRAD_H
#define TDOA_GRAD_H

#include <iostream>
#include <array>
#include "Point.h"

class grad {
public:
    // P = (x, y) аргумент
    // A B C -- неизвестные точки
    // D E F -- известные точки
    // AB_AC -- это AB - AC разность хода, аналогично для других
    // F -- целевая функция, зависящая от 6 аргументов
    // F = F(Ax, Ay, Bx, By, Cx, Cy)
    static double F(const Point &A, const Point &B, const Point &C,
                    const Point &D, const Point &E, const Point &F,
                    double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF);

    // Частная производная dF/dx
    static double dF_dx(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);

    // Частная производная dF/dy
    static double dF_dy(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);

    // градиентный метод с дроблением шага
    static Point gradMethod(const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD);
};


#endif //TDOA_GRAD_H
