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
    Solver(double alpha0, double lambda, double delta, double residual)
            : m_alpha0(alpha0), m_lambda(lambda), m_delta(delta), m_residual(residual) {}

    std::array<Point, numPoints> gradMethod(const Problem& problem);

private:
    // Начальный шаг \in (0, 1)
    const double m_alpha0;

    // Коэффициент дробления, \lambda \in (0, 1)
    const double m_lambda;

    // коэффициент, определяющий шаг \in (0, 1)
    const double m_delta;

    // Невязка
    const double m_residual;
};


#endif //TDOA_SOLVER_H
