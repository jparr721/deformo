#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {
	// Dense Vector Types
	using Vector6d = Eigen::Matrix<double, 6, 1>;
	using Vector12d = Eigen::Matrix<double, 12, 1>;

	typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorVectorXf;
	typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorVectorXi;

	// Dense Matrix Types
	using Matrix66d = Eigen::Matrix<double, 6, 6>;
	using Matrix36d = Eigen::Matrix<double, 3, 6>;
	using Matrix12d = Eigen::Matrix<double, 12, 12>;

	// Sparse Matrix Types
	using SparseMatrixXd = Eigen::SparseMatrix<double>;
} // namespace Eigen
