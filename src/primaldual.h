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
    // xが負の値をとる場合にも対応させる
    // 収束判定を入れる

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
        std::vector<T> x(this->nx, this->s0), s(this->ni, this->s0),
            ye(this->ne, T()), yi(this->ni, T()), zx(this->nx, this->s0),
            zs(this->ni, this->s0);
        std::vector<T> rpe(this->ne), rpi(this->ni), rdx(this->nx),
            rds(this->ni), rcx(this->nx), rcs(this->ni);
        std::vector<std::vector<T>> B(this->ne + this->ni,
                                      std::vector<T>(this->ne + this->ni));
        std::vector<T> ss(this->ne + this->ni);
        std::vector<T> dx(this->nx), ds(this->ni), dy(this->ne + this->ni),
            dye(this->ne), dyi(this->ni), dzx(this->nx), dzs(this->ni);
        T mu = (std::inner_product(x.begin(), x.end(), zx.begin(), T()) +
                std::inner_product(s.begin(), s.end(), zs.begin(), T())) /
               T(this->nx + this->ni);

        for (int itr = 0; itr < 20; ++itr) {
            std::cout << itr;
            for (int j = 0; j < this->nx; ++j) {
                std::cout << " " << x[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                std::cout << " " << s[j];
            }
            std::cout << std::endl;

            // STEP1
            // rp
            for (int i = 0; i < this->ne; ++i) {
                rpe[i] = _be[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpe[i] -= _Ae[i][j] * x[j];
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                rpi[i] = _bi[i] - s[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpi[i] -= _Ai[i][j] * x[j];
                }
            }
            // rd
            for (int j = 0; j < this->nx; ++j) {
                rdx[j] = _c[j] - zx[j];
                for (int i = 0; i < this->ne; ++i) {
                    rdx[j] -= _Ae[i][j] * ye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    rdx[j] -= _Ai[i][j] * yi[i];
                }
            }
            for (int j = 0; j < this->ni; ++j) {
                rds[j] = -yi[j] - zs[j];
            }
            // rc
            for (int j = 0; j < this->nx; ++j) {
                rcx[j] = this->sigma * mu - x[j] * zx[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                rcs[j] = this->sigma * mu - s[j] * zs[j];
            }

            // STEP2
            // B
            for (int i = 0; i < this->ne; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k] += _Ae[i][j] / zx[j] * x[j] * _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k + this->ne] +=
                            _Ae[i][j] / zx[j] * x[j] * _Ai[k][j];
                    }
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i + this->ne][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k] +=
                            _Ai[i][j] / zx[j] * x[j] * _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i + this->ne][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k + this->ne] +=
                            _Ai[i][j] / zx[j] * x[j] * _Ai[k][j];
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
                    ss[i] -= _Ae[i][j] * (rcx[j] - x[j] * rdx[j]) / zx[j];
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                ss[i + this->ne] = rpi[i] - (rcs[i] - s[i] * rds[i]) / zs[i];
                for (int j = 0; j < this->nx; ++j) {
                    ss[i + this->ne] -=
                        _Ai[i][j] * (rcx[j] - x[j] * rdx[j]) / zx[j];
                }
            }
            // dy
            dy = this->gauss(B, ss);
            for (int i = 0; i < this->ne; ++i) {
                dye[i] = dy[i];
            }
            for (int i = 0; i < this->ni; ++i) {
                dyi[i] = dy[i + this->ne];
            }
            // dz
            for (int j = 0; j < this->nx; ++j) {
                dzx[j] = rdx[j];
                for (int i = 0; i < this->ne; ++i) {
                    dzx[j] -= _Ae[i][j] * dye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    dzx[j] -= _Ai[i][j] * dyi[i];
                }
            }
            for (int j = 0; j < this->ni; ++j) {
                dzs[j] = rds[j] - dyi[j];
            }
            // dx
            for (int j = 0; j < this->nx; ++j) {
                dx[j] = (rcx[j] - x[j] * dzx[j]) / zx[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                ds[j] = (rcs[j] - s[j] * dzs[j]) / zs[j];
            }

            // STEP3
            T alphapx = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dx[j] < T() && -x[j] / dx[j] < alphapx) {
                    alphapx = -x[j] / dx[j];
                }
            }
            T alphaps = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (ds[j] < T() && -s[j] / ds[j] < alphaps) {
                    alphaps = -s[j] / ds[j];
                }
            }
            T alphadx = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dzx[j] < T() && -zx[j] / dzx[j] < alphadx) {
                    alphadx = -zx[j] / dzx[j];
                }
            }
            T alphads = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (dzs[j] < T() && -zs[j] / dzs[j] < alphads) {
                    alphads = -zs[j] / dzs[j];
                }
            }
            T alpha =
                this->tau * std::min({alphapx, alphaps, alphadx, alphads});

            // STEP4
            for (int j = 0; j < this->nx; ++j) {
                x[j] += alpha * dx[j];
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
                zx[j] += alpha * dzx[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                zs[j] += alpha * dzs[j];
            }
            mu = (std::inner_product(x.begin(), x.end(), zx.begin(), T()) +
                  std::inner_product(s.begin(), s.end(), zs.begin(), T())) /
                 T(this->nx + this->ni);
        }
        return x;
    }

   private:
    const int nx, ne, ni;
    T s0, sigma, tau;

    std::vector<T> gauss(std::vector<std::vector<T>> _A, std::vector<T> _b) {
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
};
}  // namespace PANSOPT
