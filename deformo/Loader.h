#pragma once

#include <Eigen/Dense>
#include <string>

class Loader {
 public:
  /*
  LoadTetgen loads faces and vertices from a tetgen file for the renderer
  */
  void LoadTetgen(Eigen::MatrixXd& V, Eigen::MatrixXd& F,
                  const std::string node_file, const std::string ele_file);

 private:
  void ReadTetgenVertexFile(Eigen::MatrixXd& V, std::ifstream& node_input);
  void ReadTetgenEleFile(Eigen::MatrixXd& F, std::ifstream& ele_input);
};
