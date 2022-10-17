#include "../src/primaldual.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDualTest, InequalityConstraint1) {
    PrimalDual<double> optimizer(2, 0, 3);
    std::vector<double> x = optimizer.UpdateVariables(
        {-4, -6}, {{}}, {}, {{1, 2}, {3, 0}, {0, 4}}, {8, 12, 12});
    // x1=4, x2=2
}

TEST(PrimalDualTest, EqualityAndInequality1) {
    PrimalDual<double> optimizer(2, 1, 2);
    std::vector<double> x = optimizer.UpdateVariables(
        {-6, -7}, {{1, -1}}, {15}, {{2, 1}, {-2, -5}}, {60, -30});
    // x1=25, x2=10
}

TEST(PrimalDualTest, NegativeDesignVariable) {
    PrimalDual<double> optimizer(2, 1, 3);
    std::vector<double> x = optimizer.UpdateVariables(
        {-2, 5}, {{4, -6}}, {30}, {{2, 8}, {-7, -5}, {-1, 0}}, {50, -10, 0});
    // x1=105/31, x2=-85/31
}
