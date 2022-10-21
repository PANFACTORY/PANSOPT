#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace PANSOPT {
template <class T>
class PrimalDual2 {
   public:
    PrimalDual2() = delete;
    PrimalDual2(const PrimalDual2<T>&) = delete;
    ~PrimalDual2() {}

    PrimalDual2(int _nx, int _ne, int _ni, T _s0 = T(1000), T _sigma = T(0.01),
                T _tau = T(0.99))
        : nx(_nx), ne(_ne), ni(_ni) {
        this->s0 = _s0;
        this->sigma = _sigma;
        this->tau = _tau;
    }

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

        for (int itr = 0; itr < 20; ++itr) {
            std::cout << itr << " ";
            T tmp = T();
            for (int j = 0; j < this->nx; ++j) {
                tmp += _c[j] * (x[j] + _l[j]);
            }
            std::cout << tmp;
            for (int j = 0; j < this->nx; ++j) {
                std::cout << " " << x[j] + _l[j];
            }
            for (int j = 0; j < this->ni; ++j) {
                std::cout << " " << s[j];
            }
            for (int j = 0; j < this->nx; ++j) {
                std::cout << " " << t[j];
            }
            std::cout << std::endl;

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
            dy = this->gauss(B, ss);
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
            T alphax = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dx[j] < T() && -x[j] / dx[j] < alphax) {
                    alphax = -x[j] / dx[j];
                }
            }
            T alphas = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (ds[j] < T() && -s[j] / ds[j] < alphas) {
                    alphas = -s[j] / ds[j];
                }
            }
            T alphat = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dt[j] < T() && -t[j] / dt[j] < alphat) {
                    alphat = -t[j] / dt[j];
                }
            }
            T alphazx = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dzx[j] < T() && -zx[j] / dzx[j] < alphazx) {
                    alphazx = -zx[j] / dzx[j];
                }
            }
            T alphazs = 1 / this->tau;
            for (int j = 0; j < this->ni; ++j) {
                if (dzs[j] < T() && -zs[j] / dzs[j] < alphazs) {
                    alphazs = -zs[j] / dzs[j];
                }
            }
            T alphazt = 1 / this->tau;
            for (int j = 0; j < this->nx; ++j) {
                if (dzt[j] < T() && -zt[j] / dzt[j] < alphazt) {
                    alphazt = -zt[j] / dzt[j];
                }
            }
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
