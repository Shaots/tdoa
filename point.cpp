#include "point.h"

double point::distance(point p1, point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}


void point::paint(point p) {
    std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
}


bool point::approximatelyEqual(double a, double b, double epsilon) {
    double diff = fabs(a - b);
    if (diff <= epsilon)
        return true;
    return diff <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int point::sgn(double val) {
    return (double(0) < val) - (val < double(0));
}

point::point(double x, double y) : x(x), y(y) {}


point &point::operator=(const point &src) {
    if (&src != this) {
        x = src.x;
        y = src.y;
    }
    return *this;
}


point &point::operator*(double alpha) {
    x *= alpha;
    y *= alpha;
    return *this;
}


bool point::operator==(const point &p) const {
    double eps = pow(10, -3);
    if (point::approximatelyEqual(x, p.x, eps)
        && point::approximatelyEqual(y, p.y, eps))
        return true;
    else
        return false;
}


