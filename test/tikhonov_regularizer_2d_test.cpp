#include "../src/tikhonov_regularizer_2d.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../src/mma.h"

using namespace PANSOPT;

TEST(TikhonovRegularizer2DTest, TikhonovRegularizer2DTest1) {
    const int nx = 10, ny = 5, nxy = nx * ny, nt = 20;
    double tau = 1e-2 * nx * ny, dt = 1e-1, w = 0.2, move = 0.01, sigma = 1e-1;

    auto index = [=](int i, int j) { return i + nx * j; };
    auto heaviside = [=](double phi) {
        if (-w <= phi && phi <= w) {
            return 0.5 + 15 / 16.0 * phi / w - 5 / 8.0 * pow(phi / w, 3.0) +
                   3 / 16.0 * pow(phi / w, 5.0);
        } else if (w < phi) {
            return 1.0;
        }
        return 0.0;
    };
    auto dheaviside = [=](double phi) {
        if (w <= phi && phi <= w) {
            return 15 / 16.0 * pow(1 - pow(phi / w, 2), 2);
        }
        return 0.0;
    };

    std::vector<double> phi(nxy, 0), chi(nxy), chi_target(nxy);
    std::vector<double> f(nt);

    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            chi_target[index(i, j)] =
                (0.3 * nx < i && i < 0.7 * nx ? 1 : 0) + sigma * dist(gen);
        }
    }

    TikhonovRegularizer2D<double> regularizer(nx, ny, tau, dt);
    MMA<double> optimizer(phi.size(), 0, 1.0, {}, {}, {},
                          std::vector<double>(phi.size(), -1.0),
                          std::vector<double>(phi.size(), 1.0));
    optimizer.move = move;

    for (int t = 0; t < nt; ++t) {
        for (int idx = 0; idx < nxy; ++idx) {
            chi[idx] = heaviside(phi[idx]);
        }

        f[t] = 0;
        std::vector<double> dfdphi(nxy);
        for (int idx = 0; idx < nxy; ++idx) {
            f[t] += pow(chi[idx] - chi_target[idx], 2.0);
            dfdphi[idx] =
                2 * (chi[idx] - chi_target[idx]) * dheaviside(phi[idx]);
        }

        optimizer.UpdateVariables(phi, f[t], regularizer.gradient(phi, dfdphi),
                                  {}, {});
    }

    auto mm = std::minmax_element(f.end() - 10, f.end());
    ASSERT_LE(*(mm.second) - *(mm.first), 1e-2);
}
