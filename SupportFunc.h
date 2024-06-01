#ifndef TDOA_SUPPORTFUNC_H
#define TDOA_SUPPORTFUNC_H

#include <array>
#include "Solver.h"
#include "Point.h"


double norm2Square(const std::array<Point, Solver::numPoints> &grad);


#endif //TDOA_SUPPORTFUNC_H
