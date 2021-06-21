#pragma once

#include <string>

#include "EigenTypes.h"
#include "MeshGenerator.h"

enum CutPlaneAxis {
    x_axis = 0,
    y_axis,
    z_axis,
};

class Mesh {
  public:
    constexpr static float kNoCutPlane = 1.01f;
    inline const static Eigen::Vector4f kMeshDefaultColor =
        Eigen::Vector4f(0.f, 0.f, 1.f, 1.f);
    inline const static Eigen::Vector4f kMeshDefaultColorInvisible =
        Eigen::Vector4f(0.f, 0.f, 1.f, 0.f);
    inline const static Eigen::Vector4f kMeshDefaultSelectedColor =
        Eigen::Vector4f(0.f, 1.f, 0.f, 0.f);
    float cut_plane = kNoCutPlane;

    std::string tetgen_flags = "zpq";

    Eigen::VectorXi sim_nodes;
    Eigen::VectorXi faces;
    Eigen::VectorXf positions;
    Eigen::VectorXf colors;

    Eigen::VectorXf rest_positions;

    Mesh(const std::string& ply_path, const float cut_plane);
    Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& T,
         float cut_plane = kNoCutPlane);

    void Update(const Eigen::VectorXf& displacements);
    void SetCutPlane(float cut_plane);
    void Reset();

    [[nodiscard]] int Size() const { return positions.rows(); }
    [[nodiscard]] int FacesSize() const { return faces.rows(); }
    [[nodiscard]] int SimNodesSize() const { return sim_nodes.rows(); }
    [[nodiscard]] int GetPositionAtFaceIndex(const int face_index) const;

  private:

    template <typename In, typename Out>
    void Vectorize(Out& out, const In& in) {
        out.resize(in.rows() * in.cols(), 1);
        In in_t = in.transpose();
        out = Eigen::Map<Out>(in_t.data(), in_t.rows() * in_t.cols());
    }

    void InitializeRenderableSurfaces(const Eigen::MatrixXf& V,
                                      const Eigen::MatrixXi& T);
    void ConstructMesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                       Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                       Eigen::MatrixXi& TT);
};
