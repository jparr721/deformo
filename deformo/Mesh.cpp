#define TETLIBRARY

#include "Mesh.h"

#include <igl/readPLY.h>

#include <Eigen/Core>
#include <filesystem>
#include <iostream>
#include <vector>

#include "Utils.h"

Mesh::Mesh(const std::string& ply_path, const float cut_plane = kNoCutPlane)
    : cut_plane(cut_plane) {
    assert(std::filesystem::path(ply_path).extension() == ".ply" &&
           "INVALID OBJECT FILE");
    Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::readPLY(ply_path, V, F);

    Eigen::MatrixXf TV;
    Eigen::MatrixXi TF;
    Eigen::MatrixXi TT;

    ConstructMesh(V, F, TV, TF, TT);

    Vectorize(sim_nodes, TT);
    Eigen::MatrixXi all_faces;
    utils::MatrixUnion(all_faces, TT, TF);

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

void Mesh::Tetrahedralize(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                          const std::string& flags, tetgenio& out) {
    tetgenio in;
    assert(MeshToTetgenio(V, F, in) && "FAILED TO CONVERT MESH TO TETGENIO");

    char* t_flags = new char[flags.size() + 1];
    try {
        std::strcpy(t_flags, flags.c_str());
        ::tetrahedralize(t_flags, &in, &out);
    } catch (int e) {
        std::cerr << __FUNCTION__ << ": Tetgen has crashed" << std::endl;
    }

    if (out.numberoftetrahedra == 0) {
        std::cerr << __FUNCTION__ << ": Tetgen failed to create tets"
                  << std::endl;
    }

    delete[] t_flags;
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
                         Eigen::MatrixXi& TT) {
    tetgenio out;

    // Generate tetrahedrals from PLC mesh with max size 1e-2.
    // const std::string tetgen_flags = "zpqa1e-1";
    const std::string tetgen_flags = "zpq";
    Tetrahedralize(V, F, tetgen_flags, out);
    assert(TetgenioToMesh(out, TV, TF, TT));
}

void Mesh::CalculateTetrahedralCoordinates(const Eigen::MatrixXf& V,
                                           const Eigen::MatrixXi& T) {
}

bool Mesh::MeshToTetgenio(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                          tetgenio& in) {
    std::vector<std::vector<float>> v;
    std::vector<std::vector<int>> f;

    utils::MatrixToList(v, V);
    utils::MatrixToList(f, F);

    // Make sure everything loaded properly.
    assert(v.size() > 0);
    assert(f.size() > 0);

    // All indices start from zero.
    in.firstnumber = 0;
    in.numberofpoints = v.size();
    in.pointlist = new double[in.numberofpoints * 3];

    // Configure points
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        assert(v[i].size() == 3);
        in.pointlist[i * 3 + 0] = v[i][0];
        in.pointlist[i * 3 + 1] = v[i][1];
        in.pointlist[i * 3 + 2] = v[i][2];
    }

    in.numberoffacets = f.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // Configure faces
    for (int i = 0; i < static_cast<int>(f.size()); ++i) {
        in.facetmarkerlist[i] = i;

        // Setup facet options
        tetgenio::facet* facet = &in.facetlist[i];
        facet->numberofpolygons = 1;
        facet->polygonlist = new tetgenio::polygon[facet->numberofpolygons];
        facet->numberofholes = 0;
        facet->holelist = NULL;
        tetgenio::polygon* polygon = &facet->polygonlist[0];

        // Setup polygon options
        polygon->numberofvertices = f[i].size();
        polygon->vertexlist = new int[polygon->numberofvertices];

        for (int j = 0; j < static_cast<int>(f[i].size()); ++j) {
            polygon->vertexlist[j] = f[i][j];
        }
    }

    return true;
}

bool Mesh::TetgenioToMesh(const tetgenio& out, Eigen::MatrixXf& V,
                          Eigen::MatrixXi& F, Eigen::MatrixXi& T) {
    std::vector<std::vector<float>> v;
    std::vector<std::vector<int>> f;
    std::vector<std::vector<int>> t;

    // Process Points
    // Crash if we don't have any points
    assert(out.pointlist != nullptr && "TETGENIO OUTPUT BUFFER EMPTY");
    v.resize(out.numberofpoints, std::vector<float>(3));
    for (int i = 0; i < out.numberofpoints; ++i) {
        v[i][0] = out.pointlist[i * 3 + 0];
        v[i][1] = out.pointlist[i * 3 + 1];
        v[i][2] = out.pointlist[i * 3 + 2];
    }

    // Process Tetrahedron
    // Crash if we don't have any tetrahedron
    assert(out.tetrahedronlist != nullptr &&
           "TETGENIO TET OUTPUT BUFFER EMPTY");
    assert(out.numberofcorners == kMaxNumCorners &&
           "WTF: INVALID NUMBER OF TETRAHEDRAL CORNERS");
    t.resize(out.numberoftetrahedra, std::vector<int>(out.numberofcorners));

    // Track out indexes to make it easier to determine if our vertex nodes
    // match.
    int min_index = 1e7;
    int max_index = -1e7;
    for (int row = 0; row < out.numberoftetrahedra; ++row) {
        for (int col = 0; col < out.numberofcorners; ++col) {
            // Index is the value representing the vertex node, with a stride of
            // 4 (num corners)
            int index = out.tetrahedronlist[row * out.numberofcorners + col];
            // Set the tetrahedron value
            t[row][col] = index;
            min_index = (min_index > index ? index : min_index);
            max_index = (max_index < index ? index : max_index);
        }
    }

    assert(min_index >= 0);
    assert(max_index >= 0);
    assert(max_index < static_cast<int>(v.size()));

    for (int row = 0; row < out.numberoftrifaces; ++row) {
        if (out.trifacemarkerlist && out.trifacemarkerlist[row] >= 0) {
            std::vector<int> face(kMaxFaceSize);
            for (int col = 0; col < kMaxFaceSize; ++col) {
                face[col] = out.trifacelist[row * kMaxFaceSize + col];
            }
            f.push_back(face);
        }
    }

    utils::ListToMatrix(V, v);
    utils::ListToMatrix(F, f);
    utils::ListToMatrix(T, t);

    return true;
}
