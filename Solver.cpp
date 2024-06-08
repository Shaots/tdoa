#include "Solver.h"
#include "SupportFunc.h"

Solver::Solver(const Point &D, const Point &E, const Point &F, double AD_BD, double AD_CD, double AE_BE, double AE_CE,
               double AF_BF, double AF_CF)
        : D(D), E(E), F(F), AD_BD(AD_BD), AD_CD(AD_CD), AE_BE(AE_BE), AE_CE(AE_CE), AF_BF(AF_BF), AF_CF(AF_CF) {}

std::array<Point, Solver::numPoints> Solver::gradMethod(double residual) {

    // Начальный шаг \in (0, 1)
    const double alpha0 = 0.9;

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
        alpha = alpha0;
        do {
            std::array<Point, numPoints> grad1 = gradient(ABC1[0], ABC1[1], ABC1[2]);

            ABC2[0] = Point::difference(ABC1[0], Point::multiply(grad1[0], alpha));
            ABC2[1] = Point::difference(ABC1[1], Point::multiply(grad1[1], alpha));
            ABC2[2] = Point::difference(ABC1[2], Point::multiply(grad1[2], alpha));
            flag2 = Fun(ABC2[0], ABC2[1], ABC2[2]) -
                    Fun(ABC1[0], ABC1[1], ABC1[2]) >
                    -alpha * delta * norm2Square(grad1);
            alpha *= lambda;
        } while (flag2);
        std::array<Point, numPoints> grad2 = gradient(ABC2[0], ABC2[1], ABC2[2]);
        flag1 = norm2Square(grad2) > residual * residual;
        ABC1 = ABC2;
    } while (flag1);
    return ABC2;
}
