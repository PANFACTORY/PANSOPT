#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace PANSOPT {

template <class T>
class Adam {
   public:
    Adam() = delete;
    Adam(const Adam<T>&) = delete;
    ~Adam() {}

    Adam(int _n, const std::vector<T>& _L, const std::vector<T>& _U,
         T _alpha = T(0.001), T _beta1 = T(0.9), T _beta2 = T(0.999),
         T _epsilon = T(1e-8))
        : n(_n), mt(_n, T()), vt(_n, T()) {
        this->alpha = _alpha;
        this->beta1 = _beta1;
        this->beta2 = _beta2;
        this->epsilon = _epsilon;
        this->t = 1;
        this->L = _L;
        this->U = _U;
    }

    void UpdateVariables(std::vector<T>& _xk, const std::vector<T>& _dfdx) {
        // バイアス補正用の係数（t回目の更新に対する補正）
        T bias_correction1 = T(1) - std::pow(this->beta1, t);
        T bias_correction2 = T(1) - std::pow(this->beta2, t);

        for (int i = 0; i < this->n; ++i) {
            // 1次モーメントの更新（慣性成分）
            this->mt[i] =
                this->beta1 * this->mt[i] + (T(1) - this->beta1) * _dfdx[i];
            // 2次モーメントの更新
            this->vt[i] = this->beta2 * this->vt[i] +
                          (T(1) - this->beta2) * pow(_dfdx[i], 2);

            // バイアス補正を行ったモーメント
            T m_hat = this->mt[i] / bias_correction1;
            T v_hat = this->vt[i] / bias_correction2;

            // 通常のAdam更新候補
            T update = this->alpha * m_hat / (std::sqrt(v_hat) + this->epsilon);
            T candidate = _xk[i] - update;

            // ボックス制約による射影とモーメンタムの補正
            if (candidate < this->L[i]) {
                candidate = this->L[i];
                // 下側の壁にぶつかった場合、負方向（下向き）の慣性があるなら打ち消す
                if (this->mt[i] < T()) {
                    this->mt[i] = T();
                }
            } else if (candidate > this->U[i]) {
                candidate = this->U[i];
                // 上側の壁にぶつかった場合、正方向（上向き）の慣性があるなら打ち消す
                if (this->mt[i] > T()) {
                    this->mt[i] = T();
                }
            }

            // 更新された値をパラメータに反映
            _xk[i] = candidate;

            ++this->t;
        }
    }

   private:
    const int n;
    int t;
    std::vector<T> mt, vt, L, U;
    T alpha, beta1, beta2, epsilon;
};
}  // namespace PANSOPT
