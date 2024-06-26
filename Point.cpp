#include "Point.h"

double Point::distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

Point Point::multiply(Point p, double alpha) {
    return Point{p.x * alpha, p.y * alpha};
}

Point Point::difference(Point p1, Point p2) {
    return Point{p1.x - p2.x, p1.y - p2.y};
}

std::string Point::toString(Point p) {
    return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")";
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