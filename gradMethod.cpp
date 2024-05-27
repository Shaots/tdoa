#include "grad.h"

double grad::F(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    return (point::distance(P, B) - point::distance(P, C) - AB_AC) *
           (point::distance(P, B) - point::distance(P, C) - AB_AC) +
           (point::distance(P, B) - point::distance(P, D) - AB_AD) *
           (point::distance(P, B) - point::distance(P, D) - AB_AD);
}

double grad::dF_dx(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double PB = point::distance(P, B);
    double PC = point::distance(P, C);
    double PD = point::distance(P, D);
    double dPB_dx = (P == B) ? point::sgn(P.x - B.x) : (P.x - B.x) / PB;
    double dPC_dx = (P == C) ? point::sgn(P.x - C.x) : (P.x - C.x) / PC;
    double dPD_dx = (P == D) ? point::sgn(P.x - D.x) : (P.x - D.x) / PD;
    return (PB - PC - AB_AC) * 2 * (dPB_dx - dPC_dx) +
           (PB - PD - AB_AD) * 2 * (dPB_dx - dPD_dx);
}

double grad::dF_dy(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double PB = point::distance(P, B);
    double PC = point::distance(P, C);
    double PD = point::distance(P, D);
    double dPB_dy = (P == B) ? point::sgn(P.y - B.y) : (P.y - B.y) / PB;
    double dPC_dy = (P == C) ? point::sgn(P.y - C.y) : (P.y - C.y) / PC;
    double dPD_dy = (P == D) ? point::sgn(P.y - D.y) : (P.y - D.y) / PD;
    return (PB - PC - AB_AC) * 2 * (dPB_dy - dPC_dy) +
           (PB - PD - AB_AD) * 2 * (dPB_dy - dPD_dy);
}

point grad::gradMethod(const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
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
    point A1{(B.x + C.x + D.x), (B.y + C.y + D.y)};
    point A2{0, 0};

    // Условие внешнего цикла
    bool flag1 = true;

    // Условие внутреннего цикла
    bool flag2 = true;
    do {
        alpha = alpha0;
        do {
            double dF_dx_x1 = dF_dx(A1, B, C, D, AB_AC, AB_AD);
            double dF_dy_x1 = dF_dy(A1, B, C, D, AB_AC, AB_AD);
            A2.x = A1.x - alpha * dF_dx_x1;
            A2.y = A1.y - alpha * dF_dy_x1;
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