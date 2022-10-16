/**
 * @file augmentedLagrange.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Augmented Lagrange method
 * @version 0.1
 * @date 2022-10-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <algorithm>

namespace PANSOPT {
/**
 * @brief Optimizer based on augmented Lagrange method
 *
 * @tparam T    Type pf variables
 */
template <class T>
class AugmentedLagrange {
   public:
    AugmentedLagrange() = delete;
    AugmentedLagrange(const AugmentedLagrange<T> &) = delete;
    ~AugmentedLagrange() {}

    /**
     * @brief Construct a new Augmented Lagrange object
     *
     * @param _n        number of design variables
     * @param _alpha
     * @param _move
     */
    AugmentedLagrange(int _n, T _alpha, T _move, T *_smin, T *_smax) : n(_n) {
        this->lambda = T();
        this->mu = T(1);
        this->alpha = _alpha;
        this->move = _move;
        this->smin = _smin;
        this->smax = _smax;
    }

    /**
     * @brief
     *
     * @param _s
     * @param _dfds
     * @param _g
     * @param _dgds
     */
    void UpdateDesignVariables(T *_s, const T *_dfds, T _g, const T *_dgds) {
        this->lambda -= this->mu * _g;
        this->mu *= this->alpha;
        for (int i = 0; i < this->n; ++i) {
            _s[i] = std::max(
                this->smin[i],
                std::min(this->smax[i],
                         _s[i] - this->move * (_dfds[i] +
                                               (this->mu * _g - this->lambda) *
                                                   _dgds[i])));
        }
    }

   private:
    const int n;
    T mu, lambda, alpha, move, *smin, *smax;
};
}  // namespace PANSOPT
