#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "EigenTypes.h"
#include "tetgen.h"

class Mesh {
 public:
  Eigen::VectorXi faces;
  Eigen::VectorXf positions;
  Eigen::VectorXi tetrahedrals;
  // For calculating cut-planes
  Eigen::MatrixXf barycenters;
  Eigen::VectorXf colors;

  Mesh(const std::string& ply_path);
  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F);
  Mesh(Eigen::MatrixXf& V, Eigen::MatrixXi& F, Eigen::MatrixXi& T);

  void Tetrahedralize(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                      const std::string& flags, tetgenio& out);

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
  constexpr static int kMaxFaceSize = 3;
  constexpr static int kMaxNumCorners = 4;
  constexpr static double kNoCutPlane = 2;

  template <typename In, typename Out>
  void Vectorize(Out& out, const In& in) {
    out.resize(in.rows() * in.cols(), 1);
    In in_t = in.transpose();
    out = Eigen::Map<Out>(in_t.data(), in_t.rows() * in_t.cols());
  }

  void InitializeRenderableSurfaces(const Eigen::MatrixXf& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& T);
  void ConstructMesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                     Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                     Eigen::MatrixXi& TT);

  void CalculateTetrahedraCoordinatesWithCutPlane(
      const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
      const Eigen::MatrixXi& T, double cut_plane = kNoCutPlane);

  bool MeshToTetgenio(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                      tetgenio& in);

  bool TetgenioToMesh(const tetgenio& out, Eigen::MatrixXf& V,
                      Eigen::MatrixXi& F, Eigen::MatrixXi& T);
  // void CheckWinding
};
