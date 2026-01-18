#pragma once
#include <vector>

namespace PANSOPT {
namespace utility {
template <class T>
class CSR {
   public:
    CSR() {}
    CSR(const std::vector<std::vector<std::pair<int, T>>>& A) {
        this->ROW = A.size();
        this->indptr = std::vector<int>(this->ROW + 1, T());
        for (int idx = 0; idx < this->ROW; ++idx) {
            this->indptr[idx + 1] = this->indptr[idx] + A[idx].size();
            std::vector<std::pair<int, T>> row = A[idx];
            std::sort(row.begin(), row.end());
            for (auto a : row) {
                this->indices.push_back(a.first);
                this->data.push_back(a.second);
            }
        }
    }

    std::vector<T> operator*(const std::vector<T>& v) const {
        std::vector<T> r(this->ROW, T());
        // #pragma omp parallel for
        for (int i = 0; i < this->ROW; ++i) {
            for (int j = this->indptr[i], jend = this->indptr[i + 1]; j < jend;
                 ++j) {
                r[i] += this->data[j] * v[this->indices[j]];
            }
        }
        return r;
    }

   private:
    int ROW = 0;
    std::vector<int> indptr, indices;
    std::vector<T> data;
};
}  // namespace utility
}  // namespace PANSOPT
