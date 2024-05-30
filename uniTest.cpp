#include "uniTest.h"

void uniTest::testOnePoint(const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD) {
    Point A = grad::gradMethod(B, C, D, AB_AC, AB_AD);
    Point::paint(A);
}

void uniTest::test() {
    const Point D = {-1.0, -1.0};
    const Point E = {1.0, 0.0};
    const Point F = {0.0, 1.0};

    // Тривиальный случай
    double AD_AE = 0;
    double AD_AF = 0;
    std::cout << "Point A = ";
    uniTest::testOnePoint(D, E, F, AD_AE, AD_AF);
    std::cout << "Answer (-1/6, -1/6)" << std::endl << std::endl;

    // Нетривиальный случай
    std::cout << "Point B = ";
    double BD_BE = sqrt(2) - sqrt(5);
    double BD_BF = sqrt(2) - 3;
    uniTest::testOnePoint(D, E, F, BD_BE, BD_BF);
    std::cout << "Answer (0, -2)" << std::endl << std::endl;

    std::cout << "Point C = ";
    double CD_CE = sqrt(2) - 3;
    double CD_CF = sqrt(2) - sqrt(5);
    uniTest::testOnePoint(D, E, F, CD_CE, CD_CF);
    std::cout << "Answer (-2, 0)" << std::endl;

}
