#include "../src/augmentedLagrange.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

using namespace PANSOPT;

TEST(AugmentedLagrangeTest, AugmentedLagrangeTest1) {
    // Solve min(x) x0 + x1 s.t. x0^2+x1^2=1
    double x[2] = {0, 0}, xmin[2] = {-2, -2}, xmax[2] = {2, 2}, f, pref,
           dfdx[2], g, dgdx[2];
    AugmentedLagrange<double> optimizer(2, 1.01, 0.01, xmin, xmax);
    for (int itr = 0; itr < 500; ++itr) {
        pref = f;
        f = x[0] + x[1];
        g = 1 - (x[0] * x[0] + x[1] * x[1]);
        std::cout << itr << " " << f << " " << g << " " << x[0] << " " << x[1]
                  << std::endl;
        if (itr > 0 && fabs(f - pref) < 1e-6 * fabs(pref) && g < 0) {
            break;
        }
        dfdx[0] = 1;
        dfdx[1] = 1;
        dgdx[0] = -2 * x[0];
        dgdx[1] = -2 * x[1];
        optimizer.UpdateDesignVariables(x, dfdx, g, dgdx);
    }
    ASSERT_LE(fabs(x[0] + 1 / sqrt(2)), 1e-5);
    ASSERT_LE(fabs(x[1] + 1 / sqrt(2)), 1e-5);
}
