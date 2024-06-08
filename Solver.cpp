#include "Solver.h"
#include "SupportFunc.h"

std::array<Point, Solver::numPoints> Solver::gradMethod(const Problem& problem) {
    // шаг
    double alpha;

    // Начальное приближение
    // [0] = A
    // [1] = B
    // [2] = C
    std::array<Point, numPoints> ABC1;
    std::array<Point, numPoints> ABC2;

    // Условие внешнего цикла
    bool flag1 = true;

    // Условие внутреннего цикла
    bool flag2 = true;
    do {
        alpha = m_alpha0;
        do {
            std::array<Point, numPoints> grad1 = problem.gradient(ABC1[0], ABC1[1], ABC1[2]);

            ABC2[0] = Point::difference(ABC1[0], Point::multiply(grad1[0], alpha));
            ABC2[1] = Point::difference(ABC1[1], Point::multiply(grad1[1], alpha));
            ABC2[2] = Point::difference(ABC1[2], Point::multiply(grad1[2], alpha));
            flag2 = problem.Fun(ABC2[0], ABC2[1], ABC2[2]) -
                    problem.Fun(ABC1[0], ABC1[1], ABC1[2]) >
                    -alpha * m_delta * norm2Square(grad1);
            alpha *= m_lambda;
        } while (flag2);
        std::array<Point, numPoints> grad2 = problem.gradient(ABC2[0], ABC2[1], ABC2[2]);
        flag1 = norm2Square(grad2) > m_residual * m_residual;
        ABC1 = ABC2;
    } while (flag1);
    return ABC2;
}
