#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {
	// Dense Vector Types
	using Vector6f = Eigen::Matrix<float, 6, 1>;
	using Vector12f = Eigen::Matrix<float, 12, 1>;

	// Dense Matrix Types
	using Matrix66f = Eigen::Matrix<float, 6, 6>;
	using Matrix36f = Eigen::Matrix<float, 3, 6>;
	using Matrix12f = Eigen::Matrix<float, 12, 12>;

	// Sparse Matrix Types
	using SparseMatrixXf = Eigen::SparseMatrix<float>;
} // namespace Eigen
