#pragma once

#include <Eigen/Dense>

void ComputeRayleighDamping(Eigen::MatrixXf& out,
                           const Eigen::MatrixXf& stiffness,
                           const Eigen::MatrixXf& mass, float mu, float lambda,
                           float modifier = 1.f);
