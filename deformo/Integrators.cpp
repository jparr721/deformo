#include "Integrators.h"

void integrators::ExplicitCentralDifference(
    Eigen::VectorXd& displacement, Eigen::Ref<const Eigen::VectorXd> forces,
    Eigen::Ref<const Eigen::MatrixXd> mass_triangular) {
  displacement = forces * mass_triangular.inverse();
}
