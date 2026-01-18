#pragma once
#include <vector>

namespace PANSOPT {
namespace utility {
template <class T, class M>
std::vector<T> ConjugateGradientSolver(const M& A, const std::vector<T>& b,
                                       T eps = T(1e-9), int itrmax = 10000) {
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
}  // namespace utility
}  // namespace PANSOPT
