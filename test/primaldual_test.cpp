#include "../src/primaldual.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDualTest, PrimalDualTest1) {
    PrimalDual<double> optimizer(2, 0, 3);
    std::vector<double> x = optimizer.UpdateVariables(
        {-4, -6}, {{}}, {}, {{1, 2}, {3, 0}, {0, 4}}, {8, 12, 12});
}

TEST(PrimalDualTest, PrimalDualTest2) {
    PrimalDual<double> optimizer(2, 1, 2);
    std::vector<double> x = optimizer.UpdateVariables(
        {-6, -7}, {{1, -1}}, {15}, {{2, 1}, {-2, -5}}, {60, -30});
}
