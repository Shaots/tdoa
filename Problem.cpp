#include "Problem.h"

double Problem::Fun(const Point& A, const Point& B, const Point& C) const {
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

std::array<Point, Problem::numPoints> Problem::gradient(const Point& A, const Point& B, const Point& C) const {
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
