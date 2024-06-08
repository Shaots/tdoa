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

    std::array<Point, numPoints> gradMethod(double residual);

    std::array<Point, numPoints> gradMethod(const Problem& problem);

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
