#include "../src/primaldual2.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDual2Test, InequalityConstraint1) {
    PrimalDual2<double> optimizer(2, 0, 1);
    std::vector<double> x = optimizer.UpdateVariables(
        {-1, 1}, {{}}, {}, {{1, 1}}, {2}, {-1, 0}, {1, 3});
    ASSERT_LE(fabs(x[0] - 1), 1e-9);
    ASSERT_LE(fabs(x[1] - 0), 1e-9);
}

TEST(PrimalDual2Test, InequalityConstraint2) {
    PrimalDual2<double> optimizer(2, 1, 0);
    std::vector<double> x = optimizer.UpdateVariables({1, -1}, {{1, 1}}, {2},
                                                      {{}}, {}, {0, 0}, {1, 3});
    ASSERT_LE(fabs(x[0] - 0), 1e-9);
    ASSERT_LE(fabs(x[1] - 2), 1e-9);
}
