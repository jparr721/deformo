#pragma once

#include <Eigen/Dense>
#include <string>

enum class CutPlaneAxis {
    x_axis = 0,
    y_axis,
    z_axis,
};

CutPlaneAxis StringToCutPlaneAxis(const std::string& input);

class Mesh {
  public:
    constexpr static float kNoCutPlane = 0.f;
    inline const static Eigen::Vector4f kMeshDefaultColor =
        Eigen::Vector4f(0.f, 0.f, 1.f, 1.f);
    inline const static Eigen::Vector4f kMeshDefaultColorInvisible =
        Eigen::Vector4f(0.f, 0.f, 1.f, 0.f);
    inline const static Eigen::Vector4f kMeshDefaultSelectedColor =
        Eigen::Vector4f(0.f, 1.f, 0.f, 0.f);

    CutPlaneAxis cut_plane_axis = CutPlaneAxis::x_axis;

    float cut_plane = kNoCutPlane;

    // Geometry Colors
    Eigen::VectorXf colors;

    // Geometry Positions
    Eigen::VectorXf positions;
    Eigen::VectorXf rest_positions;

    // Geometry Faces / Tetrahedra
    Eigen::VectorXi tetrahedral_elements;
    Eigen::VectorXi faces;

    Mesh(const std::string& ply_path, const std::string& tetgen_flags = "zpq");
    Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& T);

    void Update(const Eigen::VectorXf& displacements);

    void SetCutPlane(float cut_plane);
    void SetCutPlaneAxis(CutPlaneAxis axis);
    void SetTetgenFlags(const std::string& flags);

    void Reset();

    [[nodiscard]] int Size() const { return positions.rows(); }
    [[nodiscard]] int FacesSize() const { return faces.rows(); }
    [[nodiscard]] int TetrahedralElementsSize() const {
        return tetrahedral_elements.rows();
    }

    [[nodiscard]] static constexpr int PositionStride() { return 3; }
    [[nodiscard]] static constexpr int FacesStride() { return 4; }

    [[nodiscard]] int GetPositionAtFaceIndex(const int face_index) const;

  private:
    std::string current_file_path_;

    template <typename In, typename Out>
    void Vectorize(Out& out, const In& in) {
        out.resize(in.rows() * in.cols(), 1);
        In in_t = in.transpose();
        out = Eigen::Map<Out>(in_t.data(), in_t.rows() * in_t.cols());
    }

    void InitializeRenderableSurfaces(const Eigen::MatrixXf& V,
                                      const Eigen::MatrixXi& T);
    void ConstructMesh(Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                       Eigen::MatrixXi& TT, const Eigen::MatrixXf& V,
                       const Eigen::MatrixXi& F,
                       const std::string& tetgen_flags) const;
    Eigen::MatrixXi
    ConstructRenderedFacesFromTetrahedralElements(const Eigen::MatrixXi& F,
                                                  const Eigen::VectorXi& T);

    void InitializeFromTetgenFlagsAndFile(const std::string& ply_path,
                                          const std::string& tetgen_flags);
};
