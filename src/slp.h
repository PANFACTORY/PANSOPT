/**
 * @file slp.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Sequential Linear Programming
 * @version 0.1
 * @date 2022-10-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include "primaldual.h"

namespace PANSOPT {
/**
 * @brief Optimizer based on Sequential Linear Programming
 *
 * @tparam T    Type of variables
 */
template <class T>
class SLP {
   public:
    SLP() = delete;
    SLP(const SLP<T>&) = delete;
    ~SLP() {}

    /**
     * @brief Construct a new SLP object
     *
     * @param _n    Number of design variables
     * @param _m    Number of inequality constraint
     * @param _move Move limit
     * @param _L    Minimum value of solution vector
     * @param _U    Maximum value of solution vector
     */
    SLP(int _n, int _m, T _move, const std::vector<T>& _L,
        const std::vector<T>& _U)
        : n(_n), m(_m), solver(_n, 0, _m) {
        this->move = _move;
        this->L = _L;
        this->U = _U;
    }

    /**
     * @brief Update design variables
     *
     * @param _xk   Design variable vector
     * @param _dfdx Sensitivity vector of objective
     * @param _g    Value vector of inequality constraint
     * @param _dgdx Sensitivity matrix of inequality constraint
     */
    void UpdateVariables(std::vector<T>& _xk, const std::vector<T>& _dfdx,
                         const std::vector<T>& _g,
                         const std::vector<std::vector<T> >& _dgdx) {
        std::vector<T> l(this->n), u(this->n);
        for (int i = 0; i < this->n; ++i) {
            l[i] = std::max(this->L[i] - _xk[i], -this->move);
            u[i] = std::min(this->U[i] - _xk[i], this->move);
        }
        std::vector<T> g(this->m);
        for (int i = 0; i < this->m; ++i) {
            g[i] = -_g[i];
        }
        std::vector<T> dx =
            this->solver.UpdateVariables(_dfdx, {{}}, {}, _dgdx, g, l, u);
        for (int i = 0; i < this->n; ++i) {
            _xk[i] += dx[i];
        }
    }

   private:
    const int n, m;
    T move;
    std::vector<T> L, U;
    PrimalDual<T> solver;
};
}  // namespace PANSOPT
