#include "SupportFunc.h"

double norm2Square(const std::array<Point, Solver::numPoints> &grad) {
    return grad[0].getX() * grad[0].getX() + grad[0].getY() * grad[0].getY()
           + grad[1].getX() * grad[1].getX() + grad[1].getY() * grad[1].getY()
           + grad[2].getX() * grad[2].getX() + grad[2].getY() * grad[2].getY();
}