#include "grad.h"

double grad::F(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    return (point::distance(P, B) - point::distance(P, C) - AB_AC) *
           (point::distance(P, B) - point::distance(P, C) - AB_AC) +
           (point::distance(P, B) - point::distance(P, D) - AB_AD) *
           (point::distance(P, B) - point::distance(P, D) - AB_AD);
}

double grad::dF_dx(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double AB = point::distance(P, B);
    double AC = point::distance(P, C);
    double AD = point::distance(P, D);
    double dAB_dx = (P == B) ? 0 : (P.x - B.x) / AB;
    double dAC_dx = (P == C) ? 0 : (P.x - C.x) / AC;
    double dAD_dx = (P == D) ? 0 : (P.x - D.x) / AD;
    return (AB - AC - AB_AC) * 2 * (dAB_dx - dAC_dx) +
           (AB - AD - AB_AD) * 2 * (dAB_dx - dAD_dx);
}

double grad::dF_dy(const point &P, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка P != B, P != C, P != D
    double AB = point::distance(P, B);
    double AC = point::distance(P, C);
    double AD = point::distance(P, D);
    double dAB_dy = (P == B) ? 0 : (P.y - B.y) / AB;
    double dAC_dy = (P == C) ? 0 : (P.y - C.y) / AC;
    double dAD_dy = (P == D) ? 0 : (P.y - D.y) / AD;
    return (AB - AC - AB_AC) * 2 * (dAB_dy - dAC_dy) +
           (AB - AD - AB_AD) * 2 * (dAB_dy - dAD_dy);
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
    point A1{0, 0};
    point A2{0, 0};
    point::paint(A1);

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