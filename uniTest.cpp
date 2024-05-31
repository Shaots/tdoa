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

    std::array<Point, grad::numPoints> ABC = grad::gradMethod(D, E, F, AD_BD, AD_CD, AE_BE, AE_CE, AF_BF, AF_CF);
    std::cout << "A = ";
    Point::paint(ABC[0]);
    std::cout << "B = ";
    Point::paint(ABC[1]);
    std::cout << "C = ";
    Point::paint(ABC[2]);

    std::cout << "\nAnswer:" << std::endl;
    std::cout << "A = (-1, 1)" << std::endl;
    std::cout << "B = (1, -1)" << std::endl;
    std::cout << "C = (1, 1)" << std::endl;
}
