#ifndef TDOA_POINT_H
#define TDOA_POINT_H

#include <iostream>
#include <cmath>

class point {
public:
    static double distance(point p1, point p2);

    static void paint(point p);

public:
    point(double x, double y);

    point &operator=(const point &src);

public:
    double x;
    double y;
};


#endif //TDOA_POINT_H
