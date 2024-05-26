//
// Created by Lenovo on 26.05.2024.
//

#include "gradMethod.h"

double gradMethod::F(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    return (point::distance(A, B) - point::distance(A, C) - AB_AC) *
           (point::distance(A, B) - point::distance(A, C) - AB_AC) +
           (point::distance(A, B) - point::distance(A, D) - AB_AD) *
           (point::distance(A, B) - point::distance(A, D) - AB_AD);
}

double gradMethod::dF_dax(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка А != B, A != C, A != D
    return (point::distance(A, B) - point::distance(A, C) - AB_AC) * 2 *
           ((A.x - B.x) / point::distance(A, B) - (A.x - C.x) / point::distance(A, C)) +
           (point::distance(A, B) - point::distance(A, D) - AB_AD) * 2 *
           ((A.x - B.x) / point::distance(A, B) - (A.x - D.x) / point::distance(A, D));
}

double gradMethod::dF_day(const point &A, const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
    // проверка А != B, A != C, A != D
    return (point::distance(A, B) - point::distance(A, C) - AB_AC) * 2 *
           ((A.y - B.y) / point::distance(A, B) - (A.y - C.y) / point::distance(A, C)) +
           (point::distance(A, B) - point::distance(A, D) - AB_AD) * 2 *
           ((A.y - B.y) / point::distance(A, B) - (A.y - D.y) / point::distance(A, D));
}

point gradMethod::grad(const point &B, const point &C, const point &D, double AB_AC, double AB_AD) {
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

    point A1{100.0, 100.0};
    point A2{100.0, 100.0};
    point::paint(A1);
    bool flag1 = true;
    bool flag2 = true;
    do {
        alpha = alpha0;
        do {
            double dF_dax_x1 = dF_dax(A1, B, C, D, AB_AC, AB_AD);
            double dF_day_x1 = dF_day(A1, B, C, D, AB_AC, AB_AD);
            A2.x = A1.x - alpha * dF_dax_x1;
            A2.y = A1.y - alpha * dF_day_x1;
            flag2 = F(A2, B, C, D, AB_AC, AB_AD) - F(A1, B, C, D, AB_AC, AB_AD) >
                    -alpha * delta * (dF_dax_x1 * dF_dax_x1 + dF_day_x1 * dF_day_x1);
            alpha *= lambda;
        } while (flag2);
        flag1 = dF_dax(A2, B, C, D, AB_AC, AB_AD) * dF_dax(A2, B, C, D, AB_AC, AB_AD) +
                dF_day(A2, B, C, D, AB_AC, AB_AD) * dF_day(A2, B, C, D, AB_AC, AB_AD) > eps * eps;
        A1 = A2;
    } while (flag1);
    return A2;
}
