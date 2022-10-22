/**
 * @file primaldual.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Primal Dual method
 * @version 0.1
 * @date 2022-10-17
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace PANSOPT {
/**
 * @brief Optimizer based on Primal Dual method
 *
 * @tparam T    Type pf variables
 */
template <class T>
class PrimalDual {
   public:
    PrimalDual() = delete;
    PrimalDual(const PrimalDual<T>&) = delete;
    ~PrimalDual() {}

    /**
     * @brief Construct a new Primal Dual object
     *
     * @param _nx       Number of design variables
     * @param _ne       Number of equality constraint
     * @param _ni       Number of inequality constraint
     * @param _s0
     * @param _sigma    Decrease ratio
     * @param _tau
     */
    PrimalDual(int _nx, int _ne, int _ni, T _s0 = T(1000), T _sigma = T(0.01),
               T _tau = T(0.99))
        : nx(_nx), ne(_ne), ni(_ni) {
        this->s0 = _s0;
        this->sigma = _sigma;
        this->tau = _tau;
    }

    // TODO
    // TODO: 収束判定の閾値と反復回数のパラメータ化

    /**
     * @brief
     *
     * @param _c    Coefficient vector of objective
     * @param _Ae   Coefficient matrix of equality constraint
     * @param _be   Coefficient vector of equality constraint
     * @param _Ai   Coefficient matrix of inequality constraint
     * @param _bi   Coefficient vector of inequality constraint
     * @return std::vector<T>
     */
    std::vector<T> UpdateVariables(const std::vector<T>& _c,
                                   const std::vector<std::vector<T>>& _Ae,
                                   const std::vector<T>& _be,
                                   const std::vector<std::vector<T>>& _Ai,
                                   const std::vector<T>& _bi) {
        // STEP0
        std::vector<T> xp(this->nx, this->s0), xn(this->nx, this->s0),
            s(this->ni, this->s0), ye(this->ne, T()), yi(this->ni, T()),
            zxp(this->nx, this->s0), zxn(this->nx, this->s0),
            zs(this->ni, this->s0);
        std::vector<T> rpe(this->ne), rpi(this->ni), rdxp(this->nx),
            rdxn(this->nx), rds(this->ni), rcxp(this->nx), rcxn(this->nx),
            rcs(this->ni);
        std::vector<std::vector<T>> B(this->ne + this->ni,
                                      std::vector<T>(this->ne + this->ni));
        std::vector<T> ss(this->ne + this->ni);
        std::vector<T> dxp(this->nx), dxn(this->nx), ds(this->ni),
            dy(this->ne + this->ni), dye(this->ne), dyi(this->ni),
            dzxp(this->nx), dzxn(this->nx), dzs(this->ni);
        T mu = (std::inner_product(xp.begin(), xp.end(), zxp.begin(), T()) +
                std::inner_product(xn.begin(), xn.end(), zxn.begin(), T()) +
                std::inner_product(s.begin(), s.end(), zs.begin(), T())) /
               T(2 * this->nx + this->ni);

        for (int itr = 0; itr < 30; ++itr) {
            // STEP1
            // rp
            for (int i = 0; i < this->ne; ++i) {
                rpe[i] = _be[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpe[i] -= _Ae[i][j] * (xp[j] - xn[j]);
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                rpi[i] = _bi[i] - s[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpi[i] -= _Ai[i][j] * (xp[j] - xn[j]);
                }
            }
            // rd
            for (int j = 0; j < this->nx; ++j) {
                rdxp[j] = _c[j] - zxp[j];
                for (int i = 0; i < this->ne; ++i) {
                    rdxp[j] -= _Ae[i][j] * ye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    rdxp[j] -= _Ai[i][j] * yi[i];
                }
            }
            for (int j = 0; j < this->nx; ++j) {
                rdxn[j] = -_c[j] - zxn[j];
                for (int i = 0; i < this->ne; ++i) {
                    rdxn[j] += _Ae[i][j] * ye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    rdxn[j] += _Ai[i][j] * yi[i];
                }
            }
            for (int j = 0; j < this->ni; ++j) {
                rds[j] = -yi[j] - zs[j];
            }
            // rc
            for (int j = 0; j < this->nx; ++j) {
                rcxp[j] = this->sigma * mu - xp[j] * zxp[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                rcxn[j] = this->sigma * mu - xn[j] * zxn[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                rcs[j] = this->sigma * mu - s[j] * zs[j];
            }

            // Check convergence
            T norm_rpe = rpe.size() > 0
                             ? *std::max_element(rpe.begin(), rpe.end(),
                                                 PrimalDual<T>::comp)
                             : T();
            T norm_rpi = rpi.size() > 0
                             ? *std::max_element(rpi.begin(), rpi.end(),
                                                 PrimalDual<T>::comp)
                             : T();
            T norm_rdxp = rdxp.size() > 0
                              ? *std::max_element(rdxp.begin(), rdxp.end(),
                                                  PrimalDual<T>::comp)
                              : T();
            T norm_rdxn = rdxn.size() > 0
                              ? *std::max_element(rdxn.begin(), rdxn.end(),
                                                  PrimalDual<T>::comp)
                              : T();
            T norm_rds = rds.size() > 0
                             ? *std::max_element(rds.begin(), rds.end(),
                                                 PrimalDual<T>::comp)
                             : T();
            T norm_rcxp = rcxp.size() > 0
                              ? *std::max_element(rcxp.begin(), rcxp.end(),
                                                  PrimalDual<T>::comp)
                              : T();
            T norm_rcxn = rcxn.size() > 0
                              ? *std::max_element(rcxn.begin(), rcxn.end(),
                                                  PrimalDual<T>::comp)
                              : T();
            T norm_rcs = rcs.size() > 0
                             ? *std::max_element(rcs.begin(), rcs.end(),
                                                 PrimalDual<T>::comp)
                             : T();
            if (std::max({fabs(norm_rpe), fabs(norm_rpi), fabs(norm_rdxp),
                          fabs(norm_rdxn), fabs(norm_rds), fabs(norm_rcxp),
                          fabs(norm_rcxn), fabs(norm_rcs)}) < 1e-9) {
                break;
            }

            // STEP2
            // B
            for (int i = 0; i < this->ne; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k] += _Ae[i][j] *
                                   (xp[j] / zxp[j] + xn[j] / zxn[j]) *
                                   _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k + this->ne] +=
                            _Ae[i][j] * (xp[j] / zxp[j] + xn[j] / zxn[j]) *
                            _Ai[k][j];
                    }
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i + this->ne][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k] +=
                            _Ai[i][j] * (xp[j] / zxp[j] + xn[j] / zxn[j]) *
                            _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i + this->ne][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k + this->ne] +=
                            _Ai[i][j] * (xp[j] / zxp[j] + xn[j] / zxn[j]) *
                            _Ai[k][j];
                    }
                    if (i == k) {
                        B[i + this->ne][k + this->ne] += s[k] / zs[k];
                    }
                }
            }
            // s
            for (int i = 0; i < this->ne; ++i) {
                ss[i] = rpe[i];
                for (int j = 0; j < this->nx; ++j) {
                    ss[i] -= _Ae[i][j] * ((rcxp[j] - xp[j] * rdxp[j]) / zxp[j] -
                                          (rcxn[j] - xn[j] * rdxn[j]) / zxn[j]);
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                ss[i + this->ne] = rpi[i] - (rcs[i] - s[i] * rds[i]) / zs[i];
                for (int j = 0; j < this->nx; ++j) {
                    ss[i + this->ne] -=
                        _Ai[i][j] * ((rcxp[j] - xp[j] * rdxp[j]) / zxp[j] -
                                     (rcxn[j] - xn[j] * rdxn[j]) / zxn[j]);
                }
            }
            // dy
            dy = PrimalDual<T>::gauss(B, ss);
            for (int i = 0; i < this->ne; ++i) {
                dye[i] = dy[i];
            }
            for (int i = 0; i < this->ni; ++i) {
                dyi[i] = dy[i + this->ne];
            }
            // dz
            for (int j = 0; j < this->nx; ++j) {
                dzxp[j] = rdxp[j];
                for (int i = 0; i < this->ne; ++i) {
                    dzxp[j] -= _Ae[i][j] * dye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    dzxp[j] -= _Ai[i][j] * dyi[i];
                }
            }
            for (int j = 0; j < this->nx; ++j) {
                dzxn[j] = rdxn[j];
                for (int i = 0; i < this->ne; ++i) {
                    dzxn[j] += _Ae[i][j] * dye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    dzxn[j] += _Ai[i][j] * dyi[i];
                }
            }
            for (int j = 0; j < this->ni; ++j) {
                dzs[j] = rds[j] - dyi[j];
            }
            // dx
            for (int j = 0; j < this->nx; ++j) {
                dxp[j] = (rcxp[j] - xp[j] * dzxp[j]) / zxp[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                dxn[j] = (rcxn[j] - xn[j] * dzxn[j]) / zxn[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                ds[j] = (rcs[j] - s[j] * dzs[j]) / zs[j];
            }

            // STEP3
            T alphapxp = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dxp[j] < T() && -xp[j] / dxp[j] < alphapxp) {
                    alphapxp = -xp[j] / dxp[j];
                }
            }
            T alphapxn = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dxn[j] < T() && -xn[j] / dxn[j] < alphapxn) {
                    alphapxn = -xn[j] / dxn[j];
                }
            }
            T alphaps = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (ds[j] < T() && -s[j] / ds[j] < alphaps) {
                    alphaps = -s[j] / ds[j];
                }
            }
            T alphadxp = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dzxp[j] < T() && -zxp[j] / dzxp[j] < alphadxp) {
                    alphadxp = -zxp[j] / dzxp[j];
                }
            }
            T alphadxn = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dzxn[j] < T() && -zxn[j] / dzxn[j] < alphadxn) {
                    alphadxn = -zxn[j] / dzxn[j];
                }
            }
            T alphads = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (dzs[j] < T() && -zs[j] / dzs[j] < alphads) {
                    alphads = -zs[j] / dzs[j];
                }
            }
            T alpha = this->tau * std::min({alphapxp, alphapxn, alphaps,
                                            alphadxp, alphadxn, alphads});

            // STEP4
            for (int j = 0; j < this->nx; ++j) {
                xp[j] += alpha * dxp[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                xn[j] += alpha * dxn[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                s[j] += alpha * ds[j];
            }
            for (int i = 0; i < this->ne; ++i) {
                ye[i] += alpha * dye[i];
            }
            for (int i = 0; i < this->ni; ++i) {
                yi[i] += alpha * dyi[i];
            }
            for (int j = 0; j < this->nx; ++j) {
                zxp[j] += alpha * dzxp[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                zxn[j] += alpha * dzxn[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                zs[j] += alpha * dzs[j];
            }
            mu = (std::inner_product(xp.begin(), xp.end(), zxp.begin(), T()) +
                  std::inner_product(xn.begin(), xn.end(), zxn.begin(), T()) +
                  std::inner_product(s.begin(), s.end(), zs.begin(), T())) /
                 T(2 * this->nx + this->ni);
        }
        std::vector<T> x(this->nx);
        for (int i = 0; i < this->nx; ++i) {
            x[i] = xp[i] - xn[i];
        }
        return x;
    }

   private:
    const int nx, ne, ni;
    T s0, sigma, tau;

    static std::vector<T> gauss(std::vector<std::vector<T>> _A,
                                std::vector<T> _b) {
        for (int i = 0; i < _b.size() - 1; i++) {
            //----------Get pivot----------
            T pivot = fabs(_A[i][i]);
            int pivoti = i;
            for (int j = i + 1; j < _b.size(); j++) {
                if (pivot < fabs(_A[j][i])) {
                    pivot = fabs(_A[j][i]);
                    pivoti = j;
                }
            }

            //----------Exchange pivot----------
            if (pivoti != i) {
                std::swap(_b[i], _b[pivoti]);
                for (int j = i; j < _b.size(); j++) {
                    std::swap(_A[i][j], _A[pivoti][j]);
                }
            }

            //----------Forward erase----------
            for (int j = i + 1; j < _b.size(); j++) {
                for (int k = i + 1; k < _b.size(); k++) {
                    _A[j][k] -= _A[i][k] * _A[j][i] / _A[i][i];
                }
                _b[j] -= _b[i] * _A[j][i] / _A[i][i];
            }
        }

        //----------Back substitution----------
        std::vector<T> x(_b.size());
        for (int i = _b.size() - 1; i >= 0; i--) {
            x[i] = _b[i];
            for (int j = _b.size() - 1; j > i; j--) {
                x[i] -= x[j] * _A[i][j];
            }
            x[i] /= _A[i][i];
        }
        return x;
    }

    static bool comp(T _a, T _b) { return fabs(_a) < fabs(_b); }
};
}  // namespace PANSOPT
