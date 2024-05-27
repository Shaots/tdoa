#ifndef TDOA_GRADMETHOD_H
#define TDOA_GRADMETHOD_H

#include "point.h"

class gradMethod {
public:
    // AB_AC -- это AB - AC
    // AB_AD -- это AB - AD
    static double F(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    static double dF_dax(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    static double dF_day(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD);

    static point grad(const point &B, const point &C, const point &D, double AB_AC, double AB_AD);
};


#endif //TDOA_GRADMETHOD_H
