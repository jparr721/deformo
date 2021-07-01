#pragma once

#include "Numerics.h"

#include <Eigen/Dense>
#include <string>

namespace loader {
void ReadTetgenVertexFile(MatrixXr& V, const std::string& node_file);
void ReadTetgenFaceFile(Eigen::MatrixXi& F, const std::string& face_file);
void ReadTetgenEleFile(Eigen::MatrixXi& T, const std::string& ele_file);
} // namespace loader
