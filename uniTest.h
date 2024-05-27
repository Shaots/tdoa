#ifndef TDOA_UNITEST_H
#define TDOA_UNITEST_H

#include <iostream>
#include "point.h"
#include "grad.h"

class uniTest {
public:
    //
    static void testOnePoint(const point& B, const point& C, const point& D, double AB_AC, double AB_AD);

    static void test();
};


#endif //TDOA_UNITEST_H
