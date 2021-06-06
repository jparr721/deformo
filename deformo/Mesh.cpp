#define TETLIBRARY

#include "Mesh.h"

#include <igl/barycenter.h>
#include <igl/readPLY.h>

#include <Eigen/Core>
#include <filesystem>
#include <iostream>
#include <numeric>
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

    InitializeRenderableSurfaces(TV, TT);
}

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& T, float cut_plane)
    : cut_plane(cut_plane) {
    Vectorize(positions, V);
    Vectorize(faces, T);
    colors.resize(positions.rows());
    for (int i = 0; i < positions.rows(); i += 3) {
        colors.segment(i, 3) << 1.f, 0.f, 0.f;
    }
}

void Mesh::Update(const Eigen::VectorXf& positions_) {
    positions += positions_;
}

void Mesh::SetCutPlane(float cut_plane) {
    const int visible_elements = ceil(positions.size() * cut_plane);

    for (int i = 0; i < colors.size(); i += 3) {
        if (i < visible_elements) {
            colors.segment(i, 3) << 1.f, 0.f, 0.f;
        } else {
            colors.segment(i, 3) << 1.f, 1.f, 1.f;
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

int Mesh::GetPositionAtFaceIndex(const int face_index) const { return faces(face_index) * 3; }


void Mesh::InitializeRenderableSurfaces(const Eigen::MatrixXf& V,
                                        const Eigen::MatrixXi& T) {
    igl::barycenter(V, T, barycenters);

    CalculateTetrahedraCoordinatesWithCutPlane(V, T);
}

void Mesh::ConstructMesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                         Eigen::MatrixXf& TV, Eigen::MatrixXi& TF,
                         Eigen::MatrixXi& TT) {
    tetgenio out;

    // Generate tetrahedrals from PLC mesh with max size 1e-2.
    const std::string tetgen_flags = "zpqa1e-2";
    //const std::string tetgen_flags = "zpq";
    Tetrahedralize(V, F, tetgen_flags, out);
    assert(TetgenioToMesh(out, TV, TF, TT));
}

void Mesh::CalculateTetrahedraCoordinatesWithCutPlane(
    const Eigen::MatrixXf& V, const Eigen::MatrixXi& T) {
    assert(cut_plane <= kNoCutPlane && cut_plane >= 0);
    Eigen::VectorXf v =
        barycenters.col(2).array() - barycenters.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    std::vector<int> s;

    for (int i = 0; i < v.size(); ++i) {
        if (v(i) < cut_plane) {
            s.push_back(i);
        }
    }

    const int rows = s.size() * 4;
    constexpr int cols = 3;

    Eigen::MatrixXf V_placeholder(rows, cols);
    Eigen::MatrixXf colors_placeholder(rows, cols);
    Eigen::MatrixXi F_placeholder(rows, cols);

    for (int i = 0; i < s.size(); ++i) {
        V_placeholder.row(i * 4 + 0) = V.row(T(s.at(i), 0));
        V_placeholder.row(i * 4 + 1) = V.row(T(s.at(i), 1));
        V_placeholder.row(i * 4 + 2) = V.row(T(s.at(i), 2));
        V_placeholder.row(i * 4 + 3) = V.row(T(s.at(i), 3));

        F_placeholder.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
        F_placeholder.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
        F_placeholder.row(i * 4 + 2) << (i * 4) + 2, (i * 4) + 2, (i * 4) + 0;
        F_placeholder.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;

        // Color each face
        colors_placeholder.row(i * 4 + 0) << 1.f, 0.f, 0.f;
        colors_placeholder.row(i * 4 + 1) << 1.f, 0.f, 0.f;
        colors_placeholder.row(i * 4 + 2) << 1.f, 0.f, 0.f;
        colors_placeholder.row(i * 4 + 3) << 1.f, 0.f, 0.f;
    }

    Vectorize(positions, V_placeholder);
    Vectorize(faces, F_placeholder);
    Vectorize(colors, colors_placeholder);
}

bool Mesh::MeshToTetgenio(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
                          tetgenio& in) {
    std::vector<std::vector<float>> v;
    std::vector<std::vector<int>> f;

    utils::MatrixToList(V, v);
    utils::MatrixToList(F, f);

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

    utils::ListToMatrix(v, V);
    utils::ListToMatrix(f, F);
    utils::ListToMatrix(t, T);

    return true;
}
