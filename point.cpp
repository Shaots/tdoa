//
// Created by Lenovo on 26.05.2024.
//

#include "point.h"

double point::distance(point p1, point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

void point::paint(point p) {
    std::cout << "P: (" << p.x << ", " << p.y << ")" << std::endl;
}

point::point(double x, double y) : x(x), y(y) {}

point &point::operator=(const point &src) {
    if (&src != this) {
        x = src.x;
        y = src.y;
    }
    return *this;
}


