#include "uniTest.h"

void uniTest::test() {
    const point B = {0.0, 0.0};
    const point C = {1918.6593908214895, 4770.574154674195};
    const point D = {4111.911492926018, 3661.9046769939555};
    double AB_AC = -205.8;
    double AB_AD = -857.5;
    point A = grad::gradMethod(B, C, D, AB_AC, AB_AD);
    point::paint(A);
    std::cout << grad::F(A, B, C, D, AB_AC, AB_AD) << std::endl;
}
