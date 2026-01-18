#pragma once
#include <algorithm>
#include <numeric>
#include <vector>

#include "./utility/conjugate_gradient_solver.h"
#include "./utility/csr.h"

namespace PANSOPT {
template <class T>
class TikhonovRegularizer {
   public:
    TikhonovRegularizer(int nx, int ny, T tau, T dt)
        : nx(nx),
          ny(ny),
          nxy(nx * ny),
          tau(tau * nxy),
          dt(dt),
          elements((nx - 1) * (ny - 1)) {
        // Generate elements
        for (int i = 0; i < this->nx - 1; ++i) {
            for (int j = 0; j < this->ny - 1; ++j) {
                this->elements[i + (this->nx - 1) * j] = {
                    i + this->nx * j, (i + 1) + this->nx * j,
                    (i + 1) + this->nx * (j + 1), i + this->nx * (j + 1)};
            }
        }

        // Generate global matrix as CSR format
        T Ke[4][4] = {{cd0 * tau + cm0 / dt, cd1 * tau + cm1 / dt,
                       cd2 * tau + cm2 / dt, cd1 * tau + cm1 / dt},  //
                      {cd1 * tau + cm1 / dt, cd0 * tau + cm0 / dt,
                       cd1 * tau + cm1 / dt, cd2 * tau + cm2 / dt},
                      {cd2 * tau + cm2 / dt, cd1 * tau + cm1 / dt,
                       cd0 * tau + cm0 / dt, cd1 * tau + cm1 / dt},
                      {cd1 * tau + cm1 / dt, cd2 * tau + cm2 / dt,
                       cd1 * tau + cm1 / dt, cd0 * tau + cm0 / dt}};
        std::vector<std::vector<std::pair<int, T>>> tmp_K(this->nx * this->ny);
        for (auto e : this->elements) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    auto it = std::find_if(
                        tmp_K[e[i]].begin(), tmp_K[e[i]].end(),
                        [=](std::pair<int, T> p) { return p.first == e[j]; });
                    if (it == tmp_K[e[i]].end()) {
                        tmp_K[e[i]].push_back(
                            std::pair<int, T>(e[j], Ke[i][j]));
                    } else {
                        (*it).second += Ke[i][j];
                    }
                }
            }
        }
        this->K = utility::CSR<T>(tmp_K);
    }

    TikhonovRegularizer(const TikhonovRegularizer<T>&) = delete;
    ~TikhonovRegularizer() {}

    std::vector<T> gradient(const std::vector<T>& s_t,
                            const std::vector<T>& df) {
        T C = s_t.size() /
              std::accumulate(df.begin(), df.end(), T(),
                              [](T acc, T dfi) { return acc + fabs(dfi); });

        std::vector<T> F(this->nx * this->ny, T());
        for (auto e : this->elements) {
            T Se[4] = {-C * df[e[0]] + s_t[e[0]] / dt,  //
                       -C * df[e[1]] + s_t[e[1]] / dt,  //
                       -C * df[e[2]] + s_t[e[2]] / dt,  //
                       -C * df[e[3]] + s_t[e[3]] / dt};
            F[e[0]] += cm0 * Se[0] + cm1 * Se[1] + cm2 * Se[2] + cm1 * Se[3];
            F[e[1]] += cm1 * Se[0] + cm0 * Se[1] + cm1 * Se[2] + cm2 * Se[3];
            F[e[2]] += cm2 * Se[0] + cm1 * Se[1] + cm0 * Se[2] + cm1 * Se[3];
            F[e[3]] += cm1 * Se[0] + cm2 * Se[1] + cm1 * Se[2] + cm0 * Se[3];
        }
        std::vector<T> s_tp1 = utility::ConjugateGradientSolver(this->K, F);

        std::vector<T> gradient(s_t.size());
        for (int idx = 0; idx < s_t.size(); ++idx) {
            gradient[idx] = -(s_tp1[idx] - s_t[idx]) / dt;
        }
        return gradient;
    }

   private:
    const int nx, ny, nxy;
    const T tau, dt, cm0 = 4 / 36.0, cm1 = 2 / 36.0, cm2 = 1 / 36.0,
                     cd0 = 4 / 6.0, cd1 = -1 / 6.0, cd2 = -2 / 6.0;
    std::vector<std::vector<int>> elements;
    utility::CSR<T> K;
};
}  // namespace PANSOPT
