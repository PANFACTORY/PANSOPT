#include "../src/slp.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace PANSOPT;

TEST(SLPTest, SLPTest1) {
    // Solve min(x) x0 + x1 s.t. x0^2+x1^2<=1
    std::vector<double> x = {0, 0}, xmin = {-2, -2}, xmax = {2, 2}, dfdx(2),
                        g(1);
    std::vector<std::vector<double>> dgdx(1, std::vector<double>(2));
    double f, pref;
    SLP<double> optimizer(2, 1, 0.01, xmin, xmax);
    for (int itr = 0; itr < 500; ++itr) {
        pref = f;
        f = x[0] + x[1];
        g[0] = (x[0] * x[0] + x[1] * x[1]) - 1;
        // std::cout << itr << " " << f << " " << g[0] << " " << x[0] << " "
        //           << x[1] << std::endl;
        if (itr > 0 && fabs(f - pref) < 1e-6 * fabs(pref) && g[0] < 0) {
            break;
        }
        dfdx[0] = 1;
        dfdx[1] = 1;
        dgdx[0][0] = 2 * x[0];
        dgdx[0][1] = 2 * x[1];
        optimizer.UpdateVariables(x, dfdx, g, dgdx);
    }
    ASSERT_LE(fabs(x[0] + 1 / sqrt(2)), 1e-5);
    ASSERT_LE(fabs(x[1] + 1 / sqrt(2)), 1e-5);
}
