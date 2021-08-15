// Apparently windows defines its own min and max
#define NOMINMAX

#include "Mesh.h"

#include <igl/boundary_facets.h>
#include <igl/readPLY.h>

#include <Eigen/Core>
#include <filesystem>

#include "MeshGenerator.h"
#include "Utils.h"

auto StringToSliceAxis(const std::string& input) -> SliceAxis {
    assert(input == "X-Axis" || input == "Y-Axis" ||
           input == "Z-Axis" && "INVALID CUT PLANE AXIS INPUT");

    if (input == "X-Axis") {
        return SliceAxis::x_axis;
    }

    if (input == "Y-Axis") {
        return SliceAxis::y_axis;
    }

    if (input == "Z-Axis") {
        return SliceAxis::z_axis;
    }

    // Unreachable, just want to make the linter happy.
    return SliceAxis::x_axis;
}

Mesh::Mesh(const std::string& ply_path, const std::string& tetgen_flags)
    : current_file_path_(ply_path) {
    InitializeFromTetgenFlagsAndFile(ply_path, tetgen_flags);
}

Mesh::Mesh(const MatrixXr& V, const MatrixX<int>& T) {
    Vectorize(tetrahedral_elements, T);
    InitializeRenderableSurfaces(V, T);
}

void Mesh::Update(const VectorXr& displacements) {
    positions = displacements + rest_positions;
}

void Mesh::Reset() {
    positions = rest_positions;
    for (int i = 0; i < colors.rows(); i += 4) {
        colors.segment(i, 4) << kMeshDefaultColor;
    }
}

int Mesh::GetPositionAtFaceIndex(const int face_index) const {
    assert(face_index < tetrahedral_elements.size());
    return tetrahedral_elements(face_index) * 3;
}

void Mesh::InitializeRenderableSurfaces(const MatrixXr& V,
                                        const MatrixX<int>& T) {
    Vectorize(positions, V);
    rest_positions = positions;
    Vectorize(faces, T);
    colors.resize(V.rows() * 4);
    for (int i = 0; i < V.rows() + positions.rows(); i += 4) {
        colors.segment(i, 4) << kMeshDefaultColor;
    }
}

void Mesh::ConstructMesh(MatrixXr& TV, MatrixX<int>& TF, MatrixX<int>& TT,
                         const MatrixXr& V, const MatrixX<int>& F,
                         const std::string& tetgen_flags) const {
    tetgenio out;
    Tetrahedralize(out, V, F, tetgen_flags);
    assert(TetgenioToMesh(TV, TF, TT, out));
}

void Mesh::InitializeFromTetgenFlagsAndFile(const std::string& ply_path,
                                            const std::string& tetgen_flags) {
    assert(std::filesystem::path(ply_path).extension() == ".ply" &&
           "INVALID PLY FILE");
    MatrixXr V;
    MatrixX<int> F;
    igl::readPLY(ply_path, V, F);
    InitializeFromVerticesFacesAndTetgenFlags(V, F, tetgen_flags);
}

void Mesh::InitializeFromVerticesFacesAndTetgenFlags(
    const MatrixXr& V, const MatrixX<int>& F, const std::string& tetgen_flags) {
    MatrixXr TV;
    MatrixX<int> TF;
    MatrixX<int> TT;

    ConstructMesh(TV, TF, TT, V, F, tetgen_flags);
    Vectorize(tetrahedral_elements, TT);

    igl::boundary_facets(TT, TF);
    InitializeRenderableSurfaces(TV, TF);
}

// ================ SETTERS
void Mesh::SetSliceValue(Real value) {}

void Mesh::SetSliceAxis(const SliceAxis axis) { slice_axis = axis; }

void Mesh::SetTetgenFlags(const std::string& flags) {
    InitializeFromTetgenFlagsAndFile(current_file_path_, flags);
}

void Mesh::RefreshData(const MatrixXr& V, const MatrixX<int>& F,
                       const std::string& tetgen_flags) {
    InitializeFromVerticesFacesAndTetgenFlags(V, F, tetgen_flags);
    //InitializeRenderableSurfaces(V, F);
}
