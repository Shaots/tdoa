#ifndef TDOA_PROBLEM_H
#define TDOA_PROBLEM_H

#include <array>
#include "Point.h"

class Problem {
public:
    // Три искомой точки (A, B, C)
    static const int numPoints = 3;

    Problem(const Point &D, const Point &E, const Point &F,
           double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF)
            : D(D), E(E), F(F), AD_BD(AD_BD), AD_CD(AD_CD), AE_BE(AE_BE), AE_CE(AE_CE), AF_BF(AF_BF), AF_CF(AF_CF) {}

public:
    // A B C -- неизвестные точки
    // D E F -- известные точки
    // AB_AC -- это AB - AC разность хода, аналогично для других
    // Fun -- целевая функция, зависящая от 6 аргументов
    // Fun = Fun(Ax, Ay, Bx, By, Cx, Cy)
    double Fun(const Point &A, const Point &B, const Point &C) const;


    // Градиент функции Fun
    std::array<Point, numPoints> gradient(const Point &A, const Point &B, const Point &C) const;



private:
    const Point D;
    const Point E;
    const Point F;
    const double AD_BD;
    const double AD_CD;
    const double AE_BE;
    const double AE_CE;
    const double AF_BF;
    const double AF_CF;
};


#endif //TDOA_PROBLEM_H
