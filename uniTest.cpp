#include "uniTest.h"

void uniTest::test() {
    const Point D = {-1.0, -1.0};
    const Point E = {1.0, 0.0};
    const Point F = {0.0, 1.0};
    double AD_BD = 0;
    double AD_CD = 2 - 2 * sqrt(2);
    double AE_BE = sqrt(5) - 1;
    double AE_CE = sqrt(5) - 1;
    double AF_BF = 1 - sqrt(5);
    double AF_CF = 0;
    Problem problem(D, E, F, AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF);

    double alpha0 = 0.9;
    double lambda = 0.95;
    double delta = 0.8;
    double residual = 1e-5;
    Solver solver(alpha0, lambda, delta, residual);
    std::array<Point, Solver::numPoints> ABC = solver.gradMethod(problem);
    bool accept = true;
    const Point A = {-1.0, 1.0};
    const Point B = {1.0, -1.0};
    const Point C = {1.0, 1.0};
    accept = (ABC[0].getX() - A.getX()) * (ABC[0].getX() - A.getX()) + (ABC[0].getY() - A.getY()) * (ABC[0].getY() - A.getY()) < 1e-5
             && (ABC[1].getX() - B.getX()) * (ABC[1].getX() - B.getX()) + (ABC[1].getY() - B.getY()) * (ABC[1].getY() - B.getY()) < 1e-5
             && (ABC[2].getX() - C.getX()) * (ABC[2].getX() - C.getX()) + (ABC[2].getY() - C.getY()) * (ABC[2].getY() - C.getY()) < 1e-5;
    std::cout << "Test simple problem: " << (accept ? "Pass" : "Not pass") << std::endl;
}
