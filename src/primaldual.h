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
     * @param _s0       Initial value of xp, xn, s, zxp, zxn and zs
     * @param _sigma    Decrease ratio
     * @param _tau      Coefficient of step size
     * @param _itr_max  Number of maximum iteration
     * @param _eps      Convergence criteria
     */
    PrimalDual(int _nx, int _ne, int _ni, T _s0 = T(1000), T _sigma = T(0.01),
               T _tau = T(0.99), int _itr_max = 50, T _eps = T(1e-9))
        : nx(_nx), ne(_ne), ni(_ni), itr_max(_itr_max) {
        this->s0 = _s0;
        this->sigma = _sigma;
        this->tau = _tau;
        this->eps = _eps;
    }

    /**
     * @brief Update design variables
     *
     * @param _c    Coefficient vector of objective
     * @param _Ae   Coefficient matrix of equality constraint
     * @param _be   Coefficient vector of equality constraint
     * @param _Ai   Coefficient matrix of inequality constraint
     * @param _bi   Coefficient vector of inequality constraint
     * @return std::vector<T> Solution vector
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

        for (int itr = 0; itr < this->itr_max; ++itr) {
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
            T norm_rpe = PrimalDual<T>::vector_max(rpe);
            T norm_rpi = PrimalDual<T>::vector_max(rpi);
            T norm_rdxp = PrimalDual<T>::vector_max(rdxp);
            T norm_rdxn = PrimalDual<T>::vector_max(rdxn);
            T norm_rds = PrimalDual<T>::vector_max(rds);
            T norm_rcxp = PrimalDual<T>::vector_max(rcxp);
            T norm_rcxn = PrimalDual<T>::vector_max(rcxn);
            T norm_rcs = PrimalDual<T>::vector_max(rcs);
            if (std::max({fabs(norm_rpe), fabs(norm_rpi), fabs(norm_rdxp),
                          fabs(norm_rdxn), fabs(norm_rds), fabs(norm_rcxp),
                          fabs(norm_rcxn), fabs(norm_rcs)}) < this->eps) {
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
            T alphapxp = PrimalDual<T>::alpha_max(xp, dxp, 1 / this->tau);
            T alphapxn = PrimalDual<T>::alpha_max(xn, dxn, 1 / this->tau);
            T alphaps = PrimalDual<T>::alpha_max(s, ds, 1 / this->tau);
            T alphadxp = PrimalDual<T>::alpha_max(zxp, dzxp, 1 / this->tau);
            T alphadxn = PrimalDual<T>::alpha_max(zxn, dzxn, 1 / this->tau);
            T alphads = PrimalDual<T>::alpha_max(zs, dzs, 1 / this->tau);
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

    /**
     * @brief Update design variables with box constraint
     *
     * @param _c    Coefficient vector of objective
     * @param _Ae   Coefficient matrix of equality constraint
     * @param _be   Coefficient vector of equality constraint
     * @param _Ai   Coefficient matrix of inequality constraint
     * @param _bi   Coefficient vector of inequality constraint
     * @param _l    Minimum value of solution vector
     * @param _u    Maximum value of solution vector
     * @return std::vector<T> Solution vector
     */
    std::vector<T> UpdateVariables(const std::vector<T>& _c,
                                   const std::vector<std::vector<T>>& _Ae,
                                   const std::vector<T>& _be,
                                   const std::vector<std::vector<T>>& _Ai,
                                   const std::vector<T>& _bi,
                                   const std::vector<T>& _l,
                                   const std::vector<T>& _u) {
        // STEP0
        std::vector<T> x(this->nx, this->s0), s(this->ni, this->s0),
            t(this->nx, this->s0), ye(this->ne, T()), yi(this->ni, T()),
            zx(this->nx, this->s0), zs(this->ni, this->s0),
            zt(this->nx, this->s0);
        std::vector<T> rpe(this->ne), rpi(this->ni), rdx(this->nx),
            rds(this->ni), rdt(this->nx), rcx(this->nx), rcs(this->ni),
            rct(this->nx);
        std::vector<std::vector<T>> B(this->ne + this->ni,
                                      std::vector<T>(this->ne + this->ni));
        std::vector<T> ss(this->ne + this->ni), dy(this->ne + this->ni);
        std::vector<T> dx(this->nx), ds(this->ni), dt(this->nx), dye(this->ne),
            dyi(this->ni), dzx(this->nx), dzs(this->ni), dzt(this->nx);
        T mu = (std::inner_product(x.begin(), x.end(), zx.begin(), T()) +
                std::inner_product(s.begin(), s.end(), zs.begin(), T()) +
                std::inner_product(t.begin(), t.end(), zt.begin(), T())) /
               T(2 * this->nx + this->ni);

        for (int itr = 0; itr < this->itr_max; ++itr) {
            // STEP1
            // rp
            for (int i = 0; i < this->ne; ++i) {
                rpe[i] = -_be[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpe[i] += _Ae[i][j] * (x[j] + _l[j]);
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                rpi[i] = -_bi[i] + s[i];
                for (int j = 0; j < this->nx; ++j) {
                    rpi[i] += _Ai[i][j] * (x[j] + _l[j]);
                }
            }
            // rd
            for (int j = 0; j < this->nx; ++j) {
                rdx[j] = _c[j] - zx[j] + zt[j];
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
            for (int j = 0; j < this->nx; ++j) {
                rdt[j] = _u[j] - _l[j] - x[j] - t[j];
            }
            // rc
            for (int j = 0; j < this->nx; ++j) {
                rcx[j] = this->sigma * mu - x[j] * zx[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                rcs[j] = this->sigma * mu - s[j] * zs[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                rct[j] = this->sigma * mu - t[j] * zt[j];
            }

            // Check convergence
            T norm_rpe = PrimalDual<T>::vector_max(rpe);
            T norm_rpi = PrimalDual<T>::vector_max(rpi);
            T norm_rdx = PrimalDual<T>::vector_max(rdx);
            T norm_rds = PrimalDual<T>::vector_max(rds);
            T norm_rdt = PrimalDual<T>::vector_max(rdt);
            T norm_rcx = PrimalDual<T>::vector_max(rcx);
            T norm_rcs = PrimalDual<T>::vector_max(rcs);
            T norm_rct = PrimalDual<T>::vector_max(rct);
            if (std::max({fabs(norm_rpe), fabs(norm_rpi), fabs(norm_rdx),
                          fabs(norm_rds), fabs(norm_rdt), fabs(norm_rcx),
                          fabs(norm_rcs), fabs(norm_rct)}) < this->eps) {
                break;
            }

            // STEP2
            // B
            for (int i = 0; i < this->ne; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k] += _Ae[i][j] * x[j] /
                                   (zx[j] + x[j] * zt[j] / t[j]) * _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i][k + this->ne] += _Ae[i][j] * x[j] /
                                              (zx[j] + x[j] * zt[j] / t[j]) *
                                              _Ai[k][j];
                    }
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                for (int k = 0; k < this->ne; ++k) {
                    B[i + this->ne][k] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k] += _Ai[i][j] * x[j] /
                                              (zx[j] + x[j] * zt[j] / t[j]) *
                                              _Ae[k][j];
                    }
                }
                for (int k = 0; k < this->ni; ++k) {
                    B[i + this->ne][k + this->ne] = T();
                    for (int j = 0; j < this->nx; ++j) {
                        B[i + this->ne][k + this->ne] +=
                            _Ai[i][j] * x[j] / (zx[j] + x[j] * zt[j] / t[j]) *
                            _Ai[k][j];
                    }
                    if (i == k) {
                        B[i + this->ne][k + this->ne] += s[k] / zs[k];
                    }
                }
            }
            // s
            for (int i = 0; i < this->ne; ++i) {
                ss[i] = -rpe[i];
                for (int j = 0; j < this->nx; ++j) {
                    ss[i] -=
                        _Ae[i][j] *
                        (rcx[j] -
                         x[j] * (rdx[j] + (rct[j] - zt[j] * rdt[j]) / t[j])) /
                        (zx[j] + x[j] * zt[j] / t[j]);
                }
            }
            for (int i = 0; i < this->ni; ++i) {
                ss[i + this->ne] = -rpi[i] - (rcs[i] - s[i] * rds[i]) / zs[i];
                for (int j = 0; j < this->nx; ++j) {
                    ss[i + this->ne] -=
                        _Ai[i][j] *
                        (rcx[j] -
                         x[j] * (rdx[j] + (rct[j] - zt[j] * rdt[j]) / t[j])) /
                        (zx[j] + x[j] * zt[j] / t[j]);
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
            // dx
            for (int j = 0; j < this->nx; ++j) {
                dx[j] =
                    rcx[j] - x[j] * (rdx[j] + (rct[j] - zt[j] * rdt[j]) / t[j]);
                for (int i = 0; i < this->ne; ++i) {
                    dx[j] += x[j] * _Ae[i][j] * dye[i];
                }
                for (int i = 0; i < this->ni; ++i) {
                    dx[j] += x[j] * _Ai[i][j] * dyi[i];
                }
                dx[j] /= (zx[j] + x[j] * zt[j] / t[j]);
            }
            // ds
            for (int i = 0; i < this->ni; ++i) {
                ds[i] = (rcs[i] - s[i] * (rds[i] - dyi[i])) / zs[i];
            }
            // dt
            for (int j = 0; j < this->nx; ++j) {
                dt[j] = rdt[j] - dx[j];
            }
            // dzt
            for (int j = 0; j < this->nx; ++j) {
                dzt[j] = (rct[j] - zt[j] * dt[j]) / t[j];
            }
            // dzs
            for (int i = 0; i < this->ni; ++i) {
                dzs[i] = rds[i] - dyi[i];
            }
            // dzx
            for (int j = 0; j < this->nx; ++j) {
                dzx[j] = (rcx[j] - zx[j] * dx[j]) / x[j];
            }

            // STEP3
            T alphax = PrimalDual<T>::alpha_max(x, dx, 1 / this->tau);
            T alphas = PrimalDual<T>::alpha_max(s, ds, 1 / this->tau);
            T alphat = PrimalDual<T>::alpha_max(t, dt, 1 / this->tau);
            T alphazx = PrimalDual<T>::alpha_max(zx, dzx, 1 / this->tau);
            T alphazs = PrimalDual<T>::alpha_max(zs, dzs, 1 / this->tau);
            T alphazt = PrimalDual<T>::alpha_max(zt, dzt, 1 / this->tau);
            T alpha = this->tau * std::min({alphax, alphas, alphat, alphazx,
                                            alphazs, alphazt});

            // STEP4
            for (int j = 0; j < this->nx; ++j) {
                x[j] += alpha * dx[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                s[j] += alpha * ds[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                t[j] += alpha * dt[j];
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
            for (int j = 0; j < this->nx; ++j) {
                zt[j] += alpha * dzt[j];
            }
            mu = (std::inner_product(x.begin(), x.end(), zx.begin(), T()) +
                  std::inner_product(s.begin(), s.end(), zs.begin(), T()) +
                  std::inner_product(t.begin(), t.end(), zt.begin(), T())) /
                 T(2 * this->nx + this->ni);
        }
        for (int j = 0; j < this->nx; ++j) {
            x[j] += _l[j];
        }
        return x;
    }

   private:
    const int nx, ne, ni, itr_max;
    T s0, sigma, tau, eps;

    static std::vector<T> gauss(std::vector<std::vector<T>> _A,
                                std::vector<T> _b) {
        for (int i = 0; i < int(_b.size()) - 1; i++) {
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
        for (int i = int(_b.size()) - 1; i >= 0; i--) {
            x[i] = _b[i];
            for (int j = _b.size() - 1; j > i; j--) {
                x[i] -= x[j] * _A[i][j];
            }
            x[i] /= _A[i][i];
        }
        return x;
    }

    static T vector_max(const std::vector<T>& _v) {
        return _v.size() > 0 ? *std::max_element(_v.begin(), _v.end(),
                                                 [](T _a, T _b) {
                                                     return fabs(_a) < fabs(_b);
                                                 })
                             : T();
    }

    static T alpha_max(const std::vector<T>& _v, const std::vector<T>& _dv,
                       T _init) {
        T alpha = _init;
        for (int i = 0; i < _v.size(); ++i) {
            if (_dv[i] < T() && -_v[i] / _dv[i] < alpha) {
                alpha = -_v[i] / _dv[i];
            }
        }
        return alpha;
    }
};
}  // namespace PANSOPT
