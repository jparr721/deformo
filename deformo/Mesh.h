#pragma once

#include "Numerics.h"
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
    Real slice_value = 0;

    inline const static Vector4r kMeshDefaultColor =
        Vector4r(0.f, 0.f, 1.f, 1.f);
    inline const static Vector4r kMeshDefaultColorInvisible =
        Vector4r(0.f, 0.f, 0.f, 0.f);
    inline const static Vector4r kMeshDefaultSelectedColor =
        Vector4r(0.f, 1.f, 0.f, 1.f);

    // Geometry Colors
    VectorXr colors;

    // Geometry Positions
    VectorXr positions;
    VectorXr rest_positions;

    // Geometry Faces / Tetrahedra
    VectorX<int> tetrahedral_elements;
    VectorX<int> faces;

    Mesh(const std::string& ply_path, const std::string& tetgen_flags = "zpq");
    Mesh(const MatrixXr& V, const MatrixX<int>& T);

    void Update(const VectorXr& displacements);

    void SetSliceValue(Real value);
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

    void InitializeRenderableSurfaces(const MatrixXr& V, const MatrixX<int>& T);
    void ConstructMesh(MatrixXr& TV, MatrixX<int>& TF, MatrixX<int>& TT,
                       const MatrixXr& V, const MatrixX<int>& F,
                       const std::string& tetgen_flags) const;
    void InitializeFromTetgenFlagsAndFile(const std::string& ply_path,
                                          const std::string& tetgen_flags);
    void
    InitializeFromVerticesFacesAndTetgenFlags(const MatrixXr& V,
                                              const MatrixX<int>& F,
                                              const std::string& tetgen_flags);
};
