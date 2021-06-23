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

Mesh::Mesh(const std::string& ply_path, const float cut_plane = kNoCutPlane)
    : cut_plane(cut_plane) {
    assert(std::filesystem::path(ply_path).extension() == ".ply" &&
           "INVALID PLY FILE");
    Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::readPLY(ply_path, V, F);

    Eigen::MatrixXf TV;
    Eigen::MatrixXi TF;
    Eigen::MatrixXi TT;

    ConstructMesh(V, F, TV, TF, TT);

    for (int row = 0; row < TT.rows(); ++row) {
        const Eigen::Vector3f position_one = TV.row(TT(row, 0));
        const Eigen::Vector3f position_two = TV.row(TT(row, 1));
        const Eigen::Vector3f position_three = TV.row(TT(row, 2));
        const Eigen::Vector3f position_four = TV.row(TT(row, 3));
        const float volume = utils::ComputeTetrahedraElementVolume(
            position_one, position_two, position_three, position_four);

        std::cout << "VOLUME: " << volume << std::endl;

        assert(volume > 0);
    }

    Vectorize(sim_nodes, TT);

    // Convert TT for the matrix union
    Eigen::MatrixXi TTF;
    TTF.resize(sim_nodes.rows() / 3, 3);

    for (int i = 0; i < TTF.rows(); ++i) {
        TTF.row(i) << sim_nodes(i * 3), sim_nodes(i * 3 + 1),
            sim_nodes(i * 3 + 2);
    }

    Eigen::MatrixXi all_faces;
    utils::MatrixUnion(all_faces, TTF, TF);

    InitializeRenderableSurfaces(TV, all_faces);
}

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& T, float cut_plane)
    : cut_plane(cut_plane) {
    Vectorize(sim_nodes, T);
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
    assert(face_index < sim_nodes.size());
    return sim_nodes(face_index) * 3;
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

void Mesh::ConstructMesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                         Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                         Eigen::MatrixXi& TT) const {
    tetgenio out;
    // const std::string tetgen_flags = "zpqa1e-1";
    Tetrahedralize(out, V, F, tetgen_flags);
    assert(TetgenioToMesh(TV, TF, TT, out));
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
void Mesh::SetCutPlaneAxis(CutPlaneAxis axis) { cut_plane_axis = axis; }
void Mesh::SetTetgenFlags(const std::string& flags) { tetgen_flags = flags; }
