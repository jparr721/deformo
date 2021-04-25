#include "Integrators.h"

void integrators::ExplicitCentralDifference(
    Eigen::VectorXd& displacement, Eigen::Ref<const Eigen::VectorXd> forces,
    Eigen::Ref<const Eigen::SparseMatrixXd> M_hat) {
  // M_hat is already the inverse here so this goes smoothly.
  displacement = forces * M_hat;
}
