#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace PANSOPT {
template <class T>
class PrimalDual {
   public:
    PrimalDual() = delete;
    PrimalDual(const PrimalDual<T>&) = delete;
    ~PrimalDual() {}

    PrimalDual(int _n, int _m, T _sigma = T(0.01), T _tau = T(0.99))
        : n(_n), m(_m) {
        this->sigma = _sigma;
        this->tau = _tau;
    }

    // TODO
    // 不等式制約が扱えるように変形
    // zの初期値もパラメータ化する
    // xが負の値をとる場合にも対応させる
    // 収束判定を入れる
    void UpdateVariables(std::vector<T>& _x, const std::vector<T>& _c,
                         const std::vector<std::vector<T>>& _A,
                         const std::vector<T>& _b) {
        // STEP0
        std::vector<T> y(this->m, T()), z(this->n, 100);
        std::vector<T> rp(this->m), rd(this->n), rc(this->n);
        std::vector<std::vector<T>> B(this->m, std::vector<T>(this->m));
        std::vector<T> s(this->m);
        std::vector<T> dx(this->n), dy(this->m), dz(this->n);
        T mu = std::inner_product(_x.begin(), _x.end(), z.begin(), T()) /
               T(this->n);

        for (int itr = 0; itr < 20; ++itr) {
            std::cout << itr;
            for (int j = 0; j < this->n; ++j) {
                std::cout << " " << _x[j];
            }
            std::cout << std::endl;

            // STEP1
            for (int i = 0; i < this->m; ++i) {
                rp[i] = _b[i];
                for (int j = 0; j < this->n; ++j) {
                    rp[i] -= _A[i][j] * _x[j];
                }
            }
            for (int j = 0; j < this->n; ++j) {
                rd[j] = _c[j] - z[j];
                for (int i = 0; i < this->m; ++i) {
                    rd[j] -= _A[i][j] * y[i];
                }
            }
            for (int j = 0; j < this->n; ++j) {
                rc[j] = this->sigma * mu - _x[j] * z[j];
            }

            // STEP2
            for (int i = 0; i < this->m; ++i) {
                for (int k = 0; k < this->m; ++k) {
                    B[i][k] = T();
                    for (int j = 0; j < this->n; ++j) {
                        B[i][k] += _A[i][j] / z[j] * _x[j] * _A[k][j];
                    }
                }
            }
            for (int i = 0; i < this->m; ++i) {
                s[i] = rp[i];
                for (int j = 0; j < this->n; ++j) {
                    s[i] -= _A[i][j] * (rc[j] - _x[j] * rd[j]) / z[j];
                }
            }
            dy = this->gauss(B, s);
            for (int j = 0; j < this->n; ++j) {
                dz[j] = rd[j];
                for (int i = 0; i < this->m; ++i) {
                    dz[j] -= _A[i][j] * dy[i];
                }
            }
            for (int j = 0; j < this->n; ++j) {
                dx[j] = (rc[j] - _x[j] * dz[j]) / z[j];
            }

            // STEP3
            T alphap = 1 / this->tau;
            for (int j = 0; j < this->n; ++j) {
                if (dx[j] < T() && -_x[j] / dx[j] < alphap) {
                    alphap = -_x[j] / dx[j];
                }
            }
            T alphad = 1 / this->tau;
            for (int j = 0; j < this->n; ++j) {
                if (dz[j] < T() && -z[j] / dz[j] < alphad) {
                    alphad = -z[j] / dz[j];
                }
            }
            T alpha = this->tau * std::min(alphap, alphad);

            // STEP4
            for (int j = 0; j < this->n; ++j) {
                _x[j] += alpha * dx[j];
            }
            for (int i = 0; i < this->m; ++i) {
                y[i] += alpha * dy[i];
            }
            for (int j = 0; j < this->n; ++j) {
                z[j] += alpha * dz[j];
            }
            mu = std::inner_product(_x.begin(), _x.end(), z.begin(), T()) /
                 T(this->n);
        }
    }

   private:
    const int n, m;
    T sigma, tau;

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
