#include "Mesh.h"

#include <igl/readPLY.h>

#include <Eigen/Core>
#include <filesystem>

#include "MeshGenerator.h"
#include "Utils.h"

CutPlaneAxis StringToCutPlaneAxis(const std::string& input) {
    assert(input == "X-Axis" || input == "Y-Axis" ||
           input == "Z-Axis" && "INVALID CUT PLANE AXIS INPUT");

    if (input == "X-Axis") {
        return CutPlaneAxis::x_axis;
    }

    if (input == "Y-Axis") {
        return CutPlaneAxis::y_axis;
    }

    if (input == "Z-Axis") {
        return CutPlaneAxis::z_axis;
    }

    // Unreachable, just want to make the linter happy.
    return CutPlaneAxis::x_axis;
}

Mesh::Mesh(const std::string& ply_path, const std::string& tetgen_flags)
    : current_file_path_(ply_path) {
    InitializeFromTetgenFlagsAndFile(ply_path, tetgen_flags);
}

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& T) {
    Vectorize(tetrahedral_elements, T);
    InitializeRenderableSurfaces(V, T);
}

void Mesh::Update(const Eigen::VectorXf& displacements) {
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

void Mesh::InitializeRenderableSurfaces(const Eigen::MatrixXf& V,
                                        const Eigen::MatrixXi& T) {
    Vectorize(positions, V);
    rest_positions = positions;
    Vectorize(faces, T);
    colors.resize(V.rows() * 4);
    for (int i = 0; i < V.rows() + positions.rows(); i += 4) {
        colors.segment(i, 4) << kMeshDefaultColor;
    }
}

void Mesh::ConstructMesh(Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                         Eigen::MatrixXi& TT, const Eigen::MatrixXf& V,
                         const Eigen::MatrixXi& F,
                         const std::string& tetgen_flags) const {
    tetgenio out;
    Tetrahedralize(out, V, F, tetgen_flags);
    assert(TetgenioToMesh(TV, TF, TT, out));
}

Eigen::MatrixXi
Mesh::ConstructRenderedFacesFromTetrahedralElements(const Eigen::MatrixXi& F,
                                                    const Eigen::VectorXi& T) {
    Eigen::MatrixXi TTF;
    TTF.resize(T.rows() / 3, 3);

    for (int i = 0; i < TTF.rows(); ++i) {
        TTF.row(i) << T(i * 3), T(i * 3 + 1), T(i * 3 + 2);
    }

    Eigen::MatrixXi all_faces;
    utils::MatrixUnion(all_faces, TTF, F);
    return all_faces;
}

void Mesh::InitializeFromTetgenFlagsAndFile(const std::string& ply_path,
                                            const std::string& tetgen_flags) {
    assert(std::filesystem::path(ply_path).extension() == ".ply" &&
           "INVALID PLY FILE");
    Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::readPLY(ply_path, V, F);

    Eigen::MatrixXf TV;
    Eigen::MatrixXi TF;
    Eigen::MatrixXi TT;

    ConstructMesh(TV, TF, TT, V, F, tetgen_flags);
    Vectorize(tetrahedral_elements, TT);

    const Eigen::MatrixXi faces =
        ConstructRenderedFacesFromTetrahedralElements(TF, tetrahedral_elements);

    // Convert TT for the matrix union
    InitializeRenderableSurfaces(TV, faces);
}

// ================ SETTERS
void Mesh::SetCutPlane(float cut_plane) {
    const int visible_elements = ceil(positions.size() * cut_plane);

    for (int i = 0; i < colors.size(); i += 4) {
        if (i < visible_elements) {
            colors.segment(i, 4) << kMeshDefaultColor;
        } else {
            colors.segment(i, 4) << kMeshDefaultColorInvisible;
        }
    }
}
void Mesh::SetCutPlaneAxis(const CutPlaneAxis axis) { cut_plane_axis = axis; }

void Mesh::SetTetgenFlags(const std::string& flags) {
    InitializeFromTetgenFlagsAndFile(current_file_path_, flags);
}
