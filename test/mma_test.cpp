#include "../src/mma.h"

#include <gtest/gtest.h>

using namespace PANSOPT;

TEST(MMATest, MMATest1) {
    double C1 = 0.0624, C2 = 1.0, f = 0, g = 0;
    std::vector<double> s(5, 5.0), dfds(5, 0), dgds(5, 0);

    MMA<double> optimizer(5, 1, 1.0, {0}, {1000}, {1.0},
                          std::vector<double>(5, 1.0),
                          std::vector<double>(5, 10.0));

    for (int k = 0; k < 10; k++) {
        //  Objective
        f = C1 * (s[0] + s[1] + s[2] + s[3] + s[4]);
        dfds[0] = C1;
        dfds[1] = C1;
        dfds[2] = C1;
        dfds[3] = C1;
        dfds[4] = C1;

        //  Constraint
        g = 61.0 * pow(s[0], -3.0) + 37.0 * pow(s[1], -3.0) +
            19.0 * pow(s[2], -3.0) + 7.0 * pow(s[3], -3.0) +
            1.0 * pow(s[4], -3.0) - C2;
        dgds[0] = -3.0 * 61.0 * pow(s[0], -4.0);
        dgds[1] = -3.0 * 37.0 * pow(s[1], -4.0);
        dgds[2] = -3.0 * 19.0 * pow(s[2], -4.0);
        dgds[3] = -3.0 * 7.0 * pow(s[3], -4.0);
        dgds[4] = -3.0 * 1.0 * pow(s[4], -4.0);

        if (optimizer.IsConvergence(f)) {
            break;
        }

        optimizer.UpdateVariables(s, f, dfds, {g}, {dgds});
    }

    ASSERT_LE(fabs(s[0] - 6.01604), 1e-5);
    ASSERT_LE(fabs(s[1] - 5.30875), 1e-5);
    ASSERT_LE(fabs(s[2] - 4.49466), 1e-5);
    ASSERT_LE(fabs(s[3] - 3.50119), 1e-5);
    ASSERT_LE(fabs(s[4] - 2.15303), 1e-5);
}

TEST(MMATest, MMATest2) {
    double f = 0, g0 = 0, g1 = 0;
    std::vector<double> s = {4.0, 3.0, 2.0}, dfds(3, 0), dgds0(3, 0),
                        dgds1(3, 0);
    MMA<double> optimizer(3, 2, 1.0, {0, 0}, {1000, 1000}, {1, 1},
                          std::vector<double>(3, 0.0),
                          std::vector<double>(3, 5.0));

    for (int k = 0; k < 10; k++) {
        //  Objective
        f = pow(s[0], 2.0) + pow(s[1], 2.0) + pow(s[2], 2.0);
        dfds[0] = 2.0 * s[0];
        dfds[1] = 2.0 * s[1];
        dfds[2] = 2.0 * s[2];

        //  Constraint1
        g0 = pow(s[0] - 5.0, 2.0) + pow(s[1] - 2.0, 2.0) +
             pow(s[2] - 1.0, 2.0) - 9.0;
        dgds0[0] = 2.0 * (s[0] - 5.0);
        dgds0[1] = 2.0 * (s[1] - 2.0);
        dgds0[2] = 2.0 * (s[2] - 1.0);

        //  Constraint2
        g1 = pow(s[0] - 3.0, 2.0) + pow(s[1] - 4.0, 2.0) +
             pow(s[2] - 3.0, 2.0) - 9.0;
        dgds1[0] = 2.0 * (s[0] - 3.0);
        dgds1[1] = 2.0 * (s[1] - 4.0);
        dgds1[2] = 2.0 * (s[2] - 3.0);

        if (optimizer.IsConvergence(f)) {
            break;
        }

        optimizer.UpdateVariables(s, f, dfds, {g0, g1}, {dgds0, dgds1});
    }

    ASSERT_LE(fabs(s[0] - 2.01755), 1e-5);
    ASSERT_LE(fabs(s[1] - 1.77980), 1e-5);
    ASSERT_LE(fabs(s[2] - 1.23776), 1e-5);
}

TEST(MMATest, MMATest3) {
    double C1 = 1.0, C2 = 0.124, f = 0, g0 = 0, g1 = 0;
    std::vector<double> s = {1.5, 0.5}, dfds(2, 0), dgds0(2, 0), dgds1(2, 0);
    MMA<double> optimizer(2, 2, 1.0, {0.0, 0.0}, {1000.0, 1000.0}, {1.0, 1.0},
                          {0.2, 0.1}, {4.0, 1.6});

    for (int k = 0; k < 20; k++) {
        //  Objective
        f = C1 * s[0] * sqrt(1.0 + pow(s[1], 2.0));
        dfds[0] = C1 * sqrt(1.0 + pow(s[1], 2.0));
        dfds[1] = C1 * s[0] * s[1] / sqrt(1.0 + pow(s[1], 2.0));

        //  Constraint
        g0 = C2 * sqrt(1.0 + pow(s[1], 2.0)) *
                 (8.0 * pow(s[0], -1.0) +
                  1.0 * pow(s[0], -1.0) * pow(s[1], -1.0)) -
             1.0;
        dgds0[0] =
            C2 * sqrt(1.0 + pow(s[1], 2.0)) *
            (-8.0 * pow(s[0], -2.0) - 1.0 * pow(s[0], -2.0) * pow(s[1], -1.0));
        dgds0[1] = C2 * s[1] / sqrt(1.0 + pow(s[1], 2.0)) *
                       (8.0 * pow(s[0], -1.0) +
                        1.0 * pow(s[0], -1.0) * pow(s[1], -1.0)) +
                   C2 * sqrt(1.0 + pow(s[1], 2.0)) *
                       (-1.0 * pow(s[0], -1.0) * pow(s[1], -2.0));
        g1 = C2 * sqrt(1.0 + pow(s[1], 2.0)) *
                 (8.0 * pow(s[0], -1.0) -
                  1.0 * pow(s[0], -1.0) * pow(s[1], -1.0)) -
             1.0;
        dgds1[0] =
            C2 * sqrt(1.0 + pow(s[1], 2.0)) *
            (-8.0 * pow(s[0], -2.0) + 1.0 * pow(s[0], -2.0) * pow(s[1], -1.0));
        dgds1[1] = C2 * s[1] / sqrt(1.0 + pow(s[1], 2.0)) *
                       (8.0 * pow(s[0], -1.0) -
                        1.0 * pow(s[0], -1.0) * pow(s[1], -1.0)) +
                   C2 * sqrt(1.0 + pow(s[1], 2.0)) *
                       (1.0 * pow(s[0], -1.0) * pow(s[1], -2.0));

        if (optimizer.IsConvergence(f)) {
            break;
        }

        optimizer.UpdateVariables(s, f, dfds, {g0, g1}, {dgds0, dgds1});
    }

    ASSERT_LE(fabs(s[0] - 1.41171), 1e-5);
    ASSERT_LE(fabs(s[1] - 0.376912), 1e-5);
}
