#include "grad.h"

double grad::F(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD) {
    return (Point::distance(P, B) - Point::distance(P, C) - AB_AC) *
           (Point::distance(P, B) - Point::distance(P, C) - AB_AC) +
           (Point::distance(P, B) - Point::distance(P, D) - AB_AD) *
           (Point::distance(P, B) - Point::distance(P, D) - AB_AD);
}

double grad::dF_dx(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double PB = Point::distance(P, B);
    double PC = Point::distance(P, C);
    double PD = Point::distance(P, D);
    double dPB_dx = (P == B) ? Point::sgn(P.getX() - B.getX()) : (P.getX() - B.getX()) / PB;
    double dPC_dx = (P == C) ? Point::sgn(P.getX() - C.getX()) : (P.getX() - C.getX()) / PC;
    double dPD_dx = (P == D) ? Point::sgn(P.getX() - D.getX()) : (P.getX() - D.getX()) / PD;
    return (PB - PC - AB_AC) * 2 * (dPB_dx - dPC_dx) +
           (PB - PD - AB_AD) * 2 * (dPB_dx - dPD_dx);
}

double grad::dF_dy(const Point &P, const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double PB = Point::distance(P, B);
    double PC = Point::distance(P, C);
    double PD = Point::distance(P, D);
    double dPB_dy = (P == B) ? Point::sgn(P.getY() - B.getY()) : (P.getY() - B.getY()) / PB;
    double dPC_dy = (P == C) ? Point::sgn(P.getY() - C.getY()) : (P.getY() - C.getY()) / PC;
    double dPD_dy = (P == D) ? Point::sgn(P.getY() - D.getY()) : (P.getY() - D.getY()) / PD;
    return (PB - PC - AB_AC) * 2 * (dPB_dy - dPC_dy) +
           (PB - PD - AB_AD) * 2 * (dPB_dy - dPD_dy);
}

Point grad::gradMethod(const Point &B, const Point &C, const Point &D, double AB_AC, double AB_AD) {
    //параметр, определяющий условие окончания вычислений eps \in (0, 1)
    const double eps = pow(10, -7);

    // Начальный шаг \in (0, 1)
    const double alpha0 = 0.9;

    // шаг
    double alpha;

    // Коэффициент дробления, \lambda \in (0, 1)
    const double lambda = 0.95;

    // коэффициент, определяющий шаг \in (0, 1)
    const double delta = 0.8;

    // Начальное приближение
    Point A1{(B.getX() + C.getX() + D.getX()), (B.getY() + C.getY() + D.getY())};
    Point A2{0, 0};

    // Условие внешнего цикла
    bool flag1 = true;

    // Условие внутреннего цикла
    bool flag2 = true;
    do {
        alpha = alpha0;
        do {
            double dF_dx_x1 = dF_dx(A1, B, C, D, AB_AC, AB_AD);
            double dF_dy_x1 = dF_dy(A1, B, C, D, AB_AC, AB_AD);
            A2.setX(A1.getX() - alpha * dF_dx_x1);
            A2.setY(A1.getY() - alpha * dF_dy_x1);
            flag2 = F(A2, B, C, D, AB_AC, AB_AD) - F(A1, B, C, D, AB_AC, AB_AD) >
                    -alpha * delta * (dF_dx_x1 * dF_dx_x1 + dF_dy_x1 * dF_dy_x1);
            alpha *= lambda;
        } while (flag2);
        flag1 = dF_dx(A2, B, C, D, AB_AC, AB_AD) * dF_dx(A2, B, C, D, AB_AC, AB_AD) +
                dF_dy(A2, B, C, D, AB_AC, AB_AD) * dF_dy(A2, B, C, D, AB_AC, AB_AD) > eps * eps;
        A1 = A2;
    } while (flag1);
    return A2;
}