#ifndef TDOA_POINT_H
#define TDOA_POINT_H

#include <iostream>
#include <cmath>

class point {
public:
    static double distance(point p1, point p2);

    static void paint(point p);

    static bool approximatelyEqual(double a, double b, double epsilon);

    static int sgn(double val);

public:
    point(double x, double y);

    point &operator=(const point &src);

    point &operator*(double alpha);

    bool operator==(const point &p) const;

public:
    double x;
    double y;
};


#endif //TDOA_POINT_H
