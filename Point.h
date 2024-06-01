#ifndef TDOA_POINT_H
#define TDOA_POINT_H

#include <iostream>
#include <cmath>

class Point {
public:
    static double distance(Point p1, Point p2);

    static Point multiply(Point p, double alpha);

    static Point difference(Point p1, Point p2);

    static std::string toString(Point p);

public:
    double getX() const;

    double getY() const;

    void setX(double x_);

    void setY(double y_);


public:
    Point();

    Point(double x, double y);

    Point &operator=(const Point &src);

private:
    double x;
    double y;
};


#endif //TDOA_POINT_H
