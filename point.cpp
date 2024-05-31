#include "Point.h"

double Point::distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

Point Point::multiply(Point p, double alpha) {
    return Point{p.getX() * alpha, p.getY() * alpha};
}

Point Point::difference(Point p1, Point p2) {
    return Point{p1.getX() - p2.getX(), p1.getY() - p2.getY()};
}

void Point::paint(Point p) {
    std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
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



