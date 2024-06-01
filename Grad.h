#ifndef TDOA_GRAD_H
#define TDOA_GRAD_H

#include <iostream>
#include <array>
#include "Point.h"

class Grad {
public:
    // A B C -- неизвестные точки
    // D E F -- известные точки
    // AB_AC -- это AB - AC разность хода, аналогично для других
    // Fun -- целевая функция, зависящая от 6 аргументов
    // Fun = Fun(Ax, Ay, Bx, By, Cx, Cy)
    static double Fun(const Point &A, const Point &B, const Point &C,
                      const Point &D, const Point &E, const Point &F,
                      double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF);

    // Три искомой точки (A, B, C)
    static const int numPoints = 3;

    // Градиент функции Fun
    static std::array<Point, numPoints> gradient(const Point &A, const Point &B, const Point &C,
                                      const Point &D, const Point &E, const Point &F,
                                      double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF);

    // градиентный метод с дроблением шага
    static std::array<Point, numPoints> gradMethod(const Point &D, const Point &E, const Point &F,
                            double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF);

    // Квадрат второй нормы для градиента
    static double norm2Square(const std::array<Point, numPoints> &grad);
};


#endif //TDOA_GRAD_H
