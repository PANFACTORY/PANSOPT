#include "../src/primaldual2.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDual2Test, InequalityConstraint1) {
    PrimalDual2<double> optimizer(2, 0, 1);
    std::vector<double> x = optimizer.UpdateVariables(
        {-1, 1}, {{}}, {}, {{1, 1}}, {2}, {-1, 0}, {1, 3});
    // x1=1, x2=0
}

TEST(PrimalDual2Test, InequalityConstraint2) {
    PrimalDual2<double> optimizer(2, 1, 0);
    std::vector<double> x = optimizer.UpdateVariables({1, -1}, {{1, 1}}, {2},
                                                      {{}}, {}, {0, 0}, {1, 3});
    // x1=0, x2=2
}
