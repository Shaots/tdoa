#ifndef TDOA_SOLVER_H
#define TDOA_SOLVER_H

#include <iostream>
#include <array>
#include "Point.h"

class Solver {

    // Static
public:
    // Три искомой точки (A, B, C)
    static const int numPoints = 3;

    // Member
public:
    Solver(const Point &D, const Point &E, const Point &F,
           double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF);

    std::array<Point, numPoints> gradMethod();

private:
    // A B C -- неизвестные точки
    // D E F -- известные точки
    // AB_AC -- это AB - AC разность хода, аналогично для других
    // Fun -- целевая функция, зависящая от 6 аргументов
    // Fun = Fun(Ax, Ay, Bx, By, Cx, Cy)
    double Fun(const Point &A, const Point &B, const Point &C);


    // Градиент функции Fun
    std::array<Point, numPoints> gradient(const Point &A, const Point &B, const Point &C);

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


#endif //TDOA_SOLVER_H
