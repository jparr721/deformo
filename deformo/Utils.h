#pragma once

#include <Eigen/Dense>
#include <vector>

namespace utils {
template <typename T>
void MatrixToList(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M,
                  std::vector<std::vector<T>>& V) {
  V.resize(M.rows(), std::vector<T>(M.cols()));

  for (int row = 0; row < M.rows(); ++row) {
    for (int col = 0; col < M.cols(); ++col) {
      V[row][col] = M(row, col);
    }
  }
}

template <typename T>
void ListToMatrix(const std::vector<std::vector<T>>& V,
                  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M) {
  M.resize(V.size(), V[0].size());

  for (int row = 0; row < V.size(); ++row) {
    for (int col = 0; col < V[0].size(); ++col) {
      M(row, col) = V[row][col];
    }
  }
}
}  // namespace utils
