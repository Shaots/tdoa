#ifndef TDOA_UNITEST_H
#define TDOA_UNITEST_H

#include <iostream>
#include "Point.h"
#include "grad.h"

class uniTest {
public:
    //
    static void testOnePoint(const Point& B, const Point& C, const Point& D, double AB_AC, double AB_AD);

    static void test();
};


#endif //TDOA_UNITEST_H
