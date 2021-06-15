#pragma once

#include <string>

#include "EigenTypes.h"
#include "tetgen.h"

class Mesh {
  public:
    constexpr static float kNoCutPlane = 1.01f;
    const Eigen::Vector3f kMeshDefaultColor = Eigen::Vector3f(0.f, 0.f, 1.f);
    float cut_plane = kNoCutPlane;

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
    void Tetrahedralize(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                        const std::string& flags, tetgenio& out);
    void Reset();

    [[nodiscard]] int GetPositionAtFaceIndex(const int face_index) const;

    [[nodiscard]] int Size() const { return positions.rows(); }
    [[nodiscard]] int FacesSize() const { return faces.rows(); }
    [[nodiscard]] int SimNodesSize() const { return sim_nodes.rows(); }

  private:
    constexpr static int kMaxFaceSize = 3;
    constexpr static int kMaxNumCorners = 4;

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

    void CalculateTetrahedralCoordinates(const Eigen::MatrixXf& V,
                                         const Eigen::MatrixXi& T);

    bool MeshToTetgenio(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                        tetgenio& in);

    bool TetgenioToMesh(const tetgenio& out, Eigen::MatrixXf& V,
                        Eigen::MatrixXi& F, Eigen::MatrixXi& T);
};
