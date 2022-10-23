#pragma once
#include "primaldual.h"

namespace PANSOPT {
template <class T>
class SLP {
   public:
    SLP() = delete;
    SLP(const SLP<T>&) = delete;
    ~SLP() {}

    SLP(int _n, int _m, T _move, const std::vector<T>& _L,
        const std::vector<T>& _U)
        : n(_n), m(_m), solver(_n, 0, _m) {
        this->move = _move;
        this->L = _L;
        this->U = _U;
    }

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
