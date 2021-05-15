#pragma once

#include <Eigen/Dense>
#include <string>

namespace loader {
void ReadTetgenVertexFile(Eigen::MatrixXd& V, const std::string& node_file);
void ReadTetgenFaceFile(Eigen::MatrixXd& F, const std::string& face_file);
void ReadTetgenEleFile(Eigen::MatrixXd& T, const std::string& ele_file);
}  // namespace loader
