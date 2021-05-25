#include "Integrators.h"

void integrators::ExplicitCentralDifference(
    Eigen::VectorXf& displacement, Eigen::Ref<const Eigen::VectorXf> forces,
    Eigen::Ref<const Eigen::SparseMatrixXf> M_hat) {
  // M_hat is already the inverse here so this goes smoothly.
  displacement = forces * M_hat;
}
