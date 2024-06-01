#include "Solver.h"

double Solver::Fun(const Point &A, const Point &B, const Point &C,
                   const Point &D, const Point &E, const Point &F,
                   double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF) {

    return (Point::distance(A, D) - Point::distance(B, D) - AD_BD) *
           (Point::distance(A, D) - Point::distance(B, D) - AD_BD) +
           (Point::distance(A, D) - Point::distance(C, D) - AD_CD) *
           (Point::distance(A, D) - Point::distance(C, D) - AD_CD) +
           (Point::distance(A, E) - Point::distance(B, E) - AE_BE) *
           (Point::distance(A, E) - Point::distance(B, E) - AE_BE) +
           (Point::distance(A, E) - Point::distance(C, E) - AE_CE) *
           (Point::distance(A, E) - Point::distance(C, E) - AE_CE) +
           (Point::distance(A, F) - Point::distance(B, F) - AF_BF) *
           (Point::distance(A, F) - Point::distance(B, F) - AF_BF) +
           (Point::distance(A, F) - Point::distance(C, F) - AF_CF) *
           (Point::distance(A, F) - Point::distance(C, F) - AF_CF);
}


std::array<Point, Solver::numPoints>
Solver::gradient(const Point &A, const Point &B, const Point &C, const Point &D, const Point &E, const Point &F,
                 double AD_BD, double AD_CD, double AE_BE, double AE_CE, double AF_BF, double AF_CF) {
    std::array<Point, numPoints> grad;

    // Grad[0] = dF/dA
    double dFdAx =
            2 * (Point::distance(A, D) - Point::distance(B, D) - AD_BD) * (A.getX() - D.getX()) / Point::distance(A, D)
            +
            2 * (Point::distance(A, D) - Point::distance(C, D) - AD_CD) * (A.getX() - D.getX()) / Point::distance(A, D)
            +
            2 * (Point::distance(A, E) - Point::distance(B, E) - AE_BE) * (A.getX() - E.getX()) / Point::distance(A, E)
            +
            2 * (Point::distance(A, E) - Point::distance(C, E) - AE_CE) * (A.getX() - E.getX()) / Point::distance(A, E)
            +
            2 * (Point::distance(A, F) - Point::distance(B, F) - AF_BF) * (A.getX() - F.getX()) / Point::distance(A, F)
            +
            2 * (Point::distance(A, F) - Point::distance(C, F) - AF_CF) * (A.getX() - F.getX()) / Point::distance(A, F);

    double dFdAy =
            2 * (Point::distance(A, D) - Point::distance(B, D) - AD_BD) * (A.getY() - D.getY()) / Point::distance(A, D)
            +
            2 * (Point::distance(A, D) - Point::distance(C, D) - AD_CD) * (A.getY() - D.getY()) / Point::distance(A, D)
            +
            2 * (Point::distance(A, E) - Point::distance(B, E) - AE_BE) * (A.getY() - E.getY()) / Point::distance(A, E)
            +
            2 * (Point::distance(A, E) - Point::distance(C, E) - AE_CE) * (A.getY() - E.getY()) / Point::distance(A, E)
            +
            2 * (Point::distance(A, F) - Point::distance(B, F) - AF_BF) * (A.getY() - F.getY()) / Point::distance(A, F)
            +
            2 * (Point::distance(A, F) - Point::distance(C, F) - AF_CF) * (A.getY() - F.getY()) / Point::distance(A, F);
    grad[0].setX(dFdAx);
    grad[0].setY(dFdAy);

    // Grad[1] = dF/dB
    double dFdBx = 2 * (Point::distance(A, D) - Point::distance(B, D) - AD_BD) * (-1) * (B.getX() - D.getX()) /
                   Point::distance(B, D)
                   + 2 * (Point::distance(A, E) - Point::distance(B, E) - AE_BE) * (-1) * (B.getX() - E.getX()) /
                     Point::distance(B, E)
                   + 2 * (Point::distance(A, F) - Point::distance(B, F) - AF_BF) * (-1) * (B.getX() - F.getX()) /
                     Point::distance(B, F);

    double dFdBy = 2 * (Point::distance(A, D) - Point::distance(B, D) - AD_BD) * (-1) * (B.getY() - D.getY()) /
                   Point::distance(B, D)
                   + 2 * (Point::distance(A, E) - Point::distance(B, E) - AE_BE) * (-1) * (B.getY() - E.getY()) /
                     Point::distance(B, E)
                   + 2 * (Point::distance(A, F) - Point::distance(B, F) - AF_BF) * (-1) * (B.getY() - F.getY()) /
                     Point::distance(B, F);
    grad[1].setX(dFdBx);
    grad[1].setY(dFdBy);

    // Grad[2] = dF/dC
    double dFdCx = 2 * (Point::distance(A, D) - Point::distance(C, D) - AD_CD) * (-1) * (C.getX() - D.getX()) /
                   Point::distance(C, D)
                   + 2 * (Point::distance(A, E) - Point::distance(C, E) - AE_CE) * (-1) * (C.getX() - E.getX()) /
                     Point::distance(C, E)
                   + 2 * (Point::distance(A, F) - Point::distance(C, F) - AF_CF) * (-1) * (C.getX() - F.getX()) /
                     Point::distance(C, F);

    double dFdCy = 2 * (Point::distance(A, D) - Point::distance(C, D) - AD_CD) * (-1) * (C.getY() - D.getY()) /
                   Point::distance(C, D)
                   + 2 * (Point::distance(A, E) - Point::distance(C, E) - AE_CE) * (-1) * (C.getY() - E.getY()) /
                     Point::distance(C, E)
                   + 2 * (Point::distance(A, F) - Point::distance(C, F) - AF_CF) * (-1) * (C.getY() - F.getY()) /
                     Point::distance(C, F);
    grad[2].setX(dFdCx);
    grad[2].setY(dFdCy);
    return grad;
}

std::array<Point, Solver::numPoints> Solver::gradMethod(const Point &D, const Point &E, const Point &F,
                                                        double AD_BD, double AD_CD, double AE_BE, double AE_CE,
                                                        double AF_BF, double AF_CF) {
    //параметр, определяющий условие окончания вычислений eps \in (0, 1)
    const double eps = pow(10, -5);

    // Начальный шаг \in (0, 1)
    const double alpha0 = 0.9;

    // шаг
    double alpha;

    // Коэффициент дробления, \lambda \in (0, 1)
    const double lambda = 0.95;

    // коэффициент, определяющий шаг \in (0, 1)
    const double delta = 0.8;

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
            std::array<Point, numPoints> grad1 = Solver::gradient(ABC1[0], ABC1[1], ABC1[2], D, E, F,
                                                                  AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF);

            ABC2[0] = Point::difference(ABC1[0], Point::multiply(grad1[0], alpha));
            ABC2[1] = Point::difference(ABC1[1], Point::multiply(grad1[1], alpha));
            ABC2[2] = Point::difference(ABC1[2], Point::multiply(grad1[2], alpha));
            flag2 = Fun(ABC2[0], ABC2[1], ABC2[2], D, E, F, AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF) -
                    Fun(ABC1[0], ABC1[1], ABC1[2], D, E, F, AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF) >
                    -alpha * delta * Solver::norm2Square(grad1);
            alpha *= lambda;
        } while (flag2);
        std::array<Point, numPoints> grad2 = Solver::gradient(ABC2[0], ABC2[1], ABC2[2], D, E, F,
                                                              AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF);
        flag1 = Solver::norm2Square(grad2) > eps * eps;
        ABC1 = ABC2;
    } while (flag1);
    return ABC2;
}


double Solver::norm2Square(const std::array<Point, numPoints> &grad) {
    return grad[0].getX() * grad[0].getX() + grad[0].getY() * grad[0].getY()
         + grad[1].getX() * grad[1].getX() + grad[1].getY() * grad[1].getY()
         + grad[2].getX() * grad[2].getX() + grad[2].getY() * grad[2].getY();
}
