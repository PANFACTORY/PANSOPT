#include "../src/adam.h"

#include <gtest/gtest.h>

#include <vector>

using namespace PANSOPT;

TEST(AdamTest, AdamTest1) {
    std::vector<double> x = {1, 1, 1}, x_min = {-5, -5, -5}, x_max = {5, 5, 5};

    Adam<double> optimizer(3, x_min, x_max);

    for (int itr = 1; itr <= 8000; ++itr) {
        std::vector<double> dfdx = {2 * (x[0] + 1), 2 * (x[1] - 3),
                                    2 * (x[2] - 4)};
        optimizer.UpdateVariables(x, dfdx);

        double f = pow(x[0] + 1, 2) + pow(x[1] - 3, 2) + pow(x[2] - 4, 2);
        if (itr % 400 == 0) {
            std::cout << itr << " " << x[0] << " " << x[1] << " " << x[2] << " "
                      << f << std::endl;
        }
    }

    ASSERT_LE(fabs(x[0] - (-1)), 1e-5);
    ASSERT_LE(fabs(x[1] - 3), 1e-5);
    ASSERT_LE(fabs(x[2] - 4), 1e-5);
}

TEST(AdamTest, AdamTest2) {
    std::vector<double> x = {1, 1, 1}, x_min = {0, 0, 0}, x_max = {2, 2, 2};

    Adam<double> optimizer(3, x_min, x_max);

    for (int itr = 1; itr <= 1000; ++itr) {
        std::vector<double> dfdx = {2 * (x[0] + 1), 2 * (x[1] - 3),
                                    2 * (x[2] - 4)};
        optimizer.UpdateVariables(x, dfdx);

        double f = pow(x[0] + 1, 2) + pow(x[1] - 3, 2) + pow(x[2] - 4, 2);
        if (itr % 50 == 0) {
            std::cout << itr << " " << x[0] << " " << x[1] << " " << x[2] << " "
                      << f << std::endl;
        }
    }

    ASSERT_LE(fabs(x[0] - 0), 1e-5);
    ASSERT_LE(fabs(x[1] - 2), 1e-5);
    ASSERT_LE(fabs(x[2] - 2), 1e-5);
}
