#include "../src/primaldual.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(PrimalDualTest, InequalityConstraint1) {
    PrimalDual<double> optimizer(2, 0, 3);
    std::vector<double> x = optimizer.UpdateVariables(
        {-4, -6}, {{}}, {}, {{1, 2}, {3, 0}, {0, 4}}, {8, 12, 12});
    ASSERT_LE(fabs(x[0] - 4), 1e-9);
    ASSERT_LE(fabs(x[1] - 2), 1e-9);
}

TEST(PrimalDualTest, InequalityConstraint2) {
    PrimalDual<double> optimizer(2, 0, 1);
    std::vector<double> x = optimizer.UpdateVariables(
        {-1, 1}, {{}}, {}, {{1, 1}}, {2}, {-1, 0}, {1, 3});
    ASSERT_LE(fabs(x[0] - 1), 1e-9);
    ASSERT_LE(fabs(x[1] - 0), 1e-9);
}

TEST(PrimalDualTest, InequalityConstraint3) {
    PrimalDual<double> optimizer(2, 1, 0);
    std::vector<double> x = optimizer.UpdateVariables({1, -1}, {{1, 1}}, {2},
                                                      {{}}, {}, {0, 0}, {1, 3});
    ASSERT_LE(fabs(x[0] - 0), 1e-9);
    ASSERT_LE(fabs(x[1] - 2), 1e-9);
}

TEST(PrimalDualTest, EqualityAndInequality1) {
    PrimalDual<double> optimizer(2, 1, 2);
    std::vector<double> x = optimizer.UpdateVariables(
        {-6, -7}, {{1, -1}}, {15}, {{2, 1}, {-2, -5}}, {60, -30});
    ASSERT_LE(fabs(x[0] - 25), 1e-9);
    ASSERT_LE(fabs(x[1] - 10), 1e-9);
}

TEST(PrimalDualTest, NegativeDesignVariable) {
    PrimalDual<double> optimizer(2, 1, 3);
    std::vector<double> x = optimizer.UpdateVariables(
        {-2, 5}, {{4, -6}}, {30}, {{2, 8}, {-7, -5}, {-1, 0}}, {50, -10, 0});
    ASSERT_LE(fabs(x[0] - 105 / 31.0), 1e-9);
    ASSERT_LE(fabs(x[1] - -85 / 31.0), 1e-9);
}

TEST(rimalDualTest, Noconstraint1) {
    PrimalDual<double> optimizer(2, 0, 0);
    std::vector<double> x =
        optimizer.UpdateVariables({-1, -1}, {{}}, {}, {{}}, {}, {0, 0}, {1, 1});
    ASSERT_LE(fabs(x[0] - 1), 1e-9);
    ASSERT_LE(fabs(x[1] - 1), 1e-9);
}
