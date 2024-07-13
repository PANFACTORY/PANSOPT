/**
 * @file lsm.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Level Set Method
 * @version 0.1
 * @date 2023-07-16
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

namespace PANSOPT {
template <class T, class M>
std::vector<T> ConjugateGradientSolver(const M &A, const std::vector<T> &b,
                                       T eps = T(1e-5), int itrmax = 10000) {
    std::vector<T> xk(b.size(), T()), Axk = A * xk, rk(b.size(), T());
    std::transform(b.begin(), b.end(), Axk.begin(), rk.begin(),
                   [](T a, T b) { return a - b; });
    std::vector<T> pk = rk;
    T bnorm = sqrt(std::inner_product(b.begin(), b.end(), b.begin(), T()));
    T rkrk = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());

    for (int k = 0; k < itrmax; ++k) {
        std::vector<T> Apk = A * pk;
        T alpha =
            rkrk / std::inner_product(pk.begin(), pk.end(), Apk.begin(), T());
        std::transform(xk.begin(), xk.end(), pk.begin(), xk.begin(),
                       [=](T a, T b) { return a + alpha * b; });
        std::transform(rk.begin(), rk.end(), Apk.begin(), rk.begin(),
                       [=](T a, T b) { return a - alpha * b; });
        T rkp1rkp1 = std::inner_product(rk.begin(), rk.end(), rk.begin(), T());
        T beta = rkp1rkp1 / rkrk;
        std::transform(pk.begin(), pk.end(), rk.begin(), pk.begin(),
                       [=](T a, T b) { return beta * a + b; });
        rkrk = rkp1rkp1;
        if (sqrt(rkrk) < eps * bnorm) {
            return xk;
        }
    }

    return xk;
}

template <class T>
class CSR {
   public:
    CSR() {}
    CSR(const std::vector<std::vector<std::pair<int, T>>> &A) {
        this->ROW = A.size();
        this->indptr = std::vector<int>(this->ROW + 1, T());
        for (int idx = 0; idx < this->ROW; ++idx) {
            this->indptr[idx + 1] = this->indptr[idx] + A[idx].size();
            std::vector<std::pair<int, T>> row = A[idx];
            std::sort(row.begin(), row.end());
            for (auto a : row) {
                this->indices.push_back(a.first);
                this->data.push_back(a.second);
            }
        }
    }

    std::vector<T> operator*(const std::vector<T> &v) const {
        std::vector<T> r(this->ROW, T());
        // #pragma omp parallel for
        for (int i = 0; i < this->ROW; ++i) {
            for (int j = this->indptr[i], jend = this->indptr[i + 1]; j < jend;
                 ++j) {
                r[i] += this->data[j] * v[this->indices[j]];
            }
        }
        return r;
    }

   private:
    int ROW = 0;
    std::vector<int> indptr, indices;
    std::vector<T> data;
};

/**
 * @brief Optimizer based on Level Set Method
 *
 * @tparam T Type of variables
 */
template <class T>
class LSM {
   public:
    LSM(int nx, int ny, T tau, T dt)
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
        this->K = CSR<T>(tmp_K);
    }
    LSM(const LSM<T> &) = delete;
    ~LSM() {}

    void UpdateVariables(std::vector<T> &s, const std::vector<T> &df, T g,
                         const std::vector<T> &dg) {
        // Set augmented Lagrangian parameters
        T C = elements.size() /
              std::accumulate(df.begin(), df.end(), T(),
                              [](T acc, T dfi) { return acc + fabs(dfi); });
        T lambda = T();
        if (g >= T()) {
            T sum_df_dg =
                std::inner_product(df.begin(), df.end(), dg.begin(), T());
            T sum_dg_dg =
                std::inner_product(dg.begin(), dg.end(), dg.begin(), T());
            lambda =
                g == T() ? -sum_df_dg / sum_dg_dg : -2 * sum_df_dg / sum_dg_dg;
        }

        // Update level set function
        std::vector<T> F(this->nx * this->ny, T());
        for (auto e : this->elements) {
            T Se[4] = {-C * (df[e[0]] + lambda * dg[e[0]]) + s[e[0]] / dt,  //
                       -C * (df[e[1]] + lambda * dg[e[1]]) + s[e[1]] / dt,  //
                       -C * (df[e[2]] + lambda * dg[e[2]]) + s[e[2]] / dt,  //
                       -C * (df[e[3]] + lambda * dg[e[3]]) + s[e[3]] / dt};
            F[e[0]] += cm0 * Se[0] + cm1 * Se[1] + cm2 * Se[2] + cm1 * Se[3];
            F[e[1]] += cm1 * Se[0] + cm0 * Se[1] + cm1 * Se[2] + cm2 * Se[3];
            F[e[2]] += cm2 * Se[0] + cm1 * Se[1] + cm0 * Se[2] + cm1 * Se[3];
            F[e[3]] += cm1 * Se[0] + cm2 * Se[1] + cm1 * Se[2] + cm0 * Se[3];
        }
        s = ConjugateGradientSolver(this->K, F);
        std::transform(s.begin(), s.end(), s.begin(), [](T si) {
            return std::max(T(-1), std::min(T(1), si));
        });
    }

   private:
    const int nx, ny, nxy;
    const T tau, dt, cm0 = 4 / 36.0, cm1 = 2 / 36.0, cm2 = 1 / 36.0,
                     cd0 = 4 / 6.0, cd1 = -1 / 6.0, cd2 = -2 / 6.0;
    std::vector<std::vector<int>> elements;
    CSR<T> K;
};
}  // namespace PANSOPT
