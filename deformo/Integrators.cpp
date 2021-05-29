#include "Integrators.h"

void integrators::ExplicitCentralDifference(
    Eigen::VectorXf& displacement, const Eigen::VectorXf& forces,
    const Eigen::FullPivLU<Eigen::MatrixXf>& M_hat) {
  // M_hat is already in triangular form, we just need to solve
  displacement = M_hat.solve(forces);
}
