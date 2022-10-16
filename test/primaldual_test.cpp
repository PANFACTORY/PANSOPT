#include "../src/primaldual.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDualTest, PrimalDualTest1) {
    PrimalDual<double> optimizer(5, 3);
    std::vector<double> x = {100, 100, 100, 100, 100};
    std::vector<std::vector<double>> A = {
        {1, 2, 1, 0, 0}, {3, 0, 0, 1, 0}, {0, 4, 0, 0, 1}};
    std::vector<double> b = {8, 12, 12}, c = {-4, -6, 0, 0, 0};
    optimizer.UpdateVariables(x, c, A, b);
}
