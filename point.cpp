#include "Point.h"

double Point::distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}


void Point::paint(Point p) {
    std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
}


bool Point::approximatelyEqual(double a, double b, double epsilon) {
    double diff = fabs(a - b);
    if (diff <= epsilon)
        return true;
    return diff <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int Point::sgn(double val) {
    return (double(0) < val) - (val < double(0));
}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

void Point::setX(double x_) {
    x = x_;
}

void Point::setY(double y_) {
    y = y_;
}

Point::Point() {
    x = 0;
    y = 0;
}

Point::Point(double x, double y) : x(x), y(y) {}


Point &Point::operator=(const Point &src) {
    if (&src != this) {
        x = src.x;
        y = src.y;
    }
    return *this;
}


Point &Point::operator*(double alpha) {
    x *= alpha;
    y *= alpha;
    return *this;
}


bool Point::operator==(const Point &p) const {
    double eps = pow(10, -3);
    if (Point::approximatelyEqual(x, p.x, eps)
        && Point::approximatelyEqual(y, p.y, eps))
        return true;
    else
        return false;
}