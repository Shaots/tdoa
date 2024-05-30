#ifndef TDOA_POINT_H
#define TDOA_POINT_H

#include <iostream>
#include <cmath>

class Point {
public:
    static double distance(Point p1, Point p2);

    static void paint(Point p);

    static bool approximatelyEqual(double a, double b, double epsilon);

    static int sgn(double val);

public:
    double getX() const;

    double getY() const;

    void setX(double x_);

    void setY(double y_);


public:
    Point(double x, double y);

    Point &operator=(const Point &src);

    Point &operator*(double alpha);

    bool operator==(const Point &p) const;

private:
    double x;
    double y;
};


#endif //TDOA_POINT_H
