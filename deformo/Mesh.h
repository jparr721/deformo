#pragma once

#include <Eigen/Dense>
#include <string>

enum class SliceAxis {
    x_axis = 0,
    y_axis,
    z_axis,
};

auto StringToSliceAxis(const std::string& input) -> SliceAxis;

class Mesh {
  public:
    SliceAxis slice_axis = SliceAxis::x_axis;
    float slice_value = 0;

    inline const static Eigen::Vector4f kMeshDefaultColor =
        Eigen::Vector4f(0.f, 0.f, 1.f, 1.f);
    inline const static Eigen::Vector4f kMeshDefaultColorInvisible =
        Eigen::Vector4f(0.f, 0.f, 0.f, 0.f);
    inline const static Eigen::Vector4f kMeshDefaultSelectedColor =
        Eigen::Vector4f(0.f, 1.f, 0.f, 1.f);

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

    void SetSliceValue(float value);
    void SetSliceAxis(SliceAxis axis);
    void SetTetgenFlags(const std::string& flags);

    void Reset();

    [[nodiscard]] auto Size() const -> int { return positions.rows(); }
    [[nodiscard]] auto FacesSize() const -> int { return faces.rows(); }
    [[nodiscard]] auto TetrahedralElementsSize() const -> int {
        return tetrahedral_elements.rows();
    }
    [[nodiscard]] static constexpr auto PositionStride() -> int { return 3; }
    [[nodiscard]] static constexpr auto FacesStride() -> int { return 4; }

    [[nodiscard]] auto GetPositionAtFaceIndex(const int face_index) const
        -> int;

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
    void InitializeFromTetgenFlagsAndFile(const std::string& ply_path,
                                          const std::string& tetgen_flags);
    auto ConstructRenderedFacesFromTetrahedralElements(const Eigen::MatrixXi& F,
                                                       const Eigen::VectorXi& T)
        -> Eigen::MatrixXi;
};
