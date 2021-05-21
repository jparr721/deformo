#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "EigenTypes.h"
#include "tetgen.h"

class Mesh {
 public:
  Eigen::RowMajorVectorXi faces;
  Eigen::RowMajorVectorXf positions;
  Eigen::RowMajorVectorXi tetrahedrals;
  Eigen::RowMajorVectorXf colors;

  Mesh(const std::string& ply_path);
  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F);
  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
       const Eigen::MatrixXi& T);

  void Tetrahedralize(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                      const std::string& flags);

  [[nodiscard]] const float* data() { return positions.data(); }
  [[nodiscard]] int size_bytes() { return positions.size() * sizeof(float); }
  [[nodiscard]] int size() { return positions.size(); }

  [[nodiscard]] const float* colors_data() { return colors.data(); }
  [[nodiscard]] int colors_size_bytes() {
    return colors.size() * sizeof(float);
  }
  [[nodiscard]] int colors_size() { return colors.size(); }

  [[nodiscard]] const int* faces_data() { return faces.data(); }
  [[nodiscard]] int faces_size_bytes() { return faces.size() * sizeof(int); }
  [[nodiscard]] int faces_size() { return faces.size(); }

 private:
  template <typename Map, typename In, typename Out>
  void Vectorize(Out& out, const In& in) {
    out.resize(in.rows() * in.cols(), 1);
    In in_t = in.transpose();
    out = Eigen::Map<Map>(in_t.data(), in_t.rows() * in_t.cols());
  }

  void GetTetrahedralCombinations(std::vector<std::vector<int>>& combinations,
                                  const Eigen::VectorXi& tet);

  bool MeshToTetgenio(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                      tetgenio& in);

  bool TetgenioToMesh(const tetgenio& out, Eigen::MatrixXf& V,
                      Eigen::MatrixXi& F, Eigen::MatrixXi& T);
  // void CheckWinding
};
