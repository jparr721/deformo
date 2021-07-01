#define TETLIBRARY

#include "MeshGenerator.h"
#include "Utils.h"

#include <cassert>
#include <iostream>
#include <vector>

namespace {
constexpr int kMaxFaceSize = 3;
constexpr int kMaxNumCorners = 4;
} // namespace

bool TetgenioToMesh(MatrixXr& V, Eigen::MatrixXi& F, Eigen::MatrixXi& T,
                    const tetgenio& out) {
    std::vector<std::vector<Real>> v;
    std::vector<std::vector<int>> f;
    std::vector<std::vector<int>> t;

    // Process Points
    // Crash if we don't have any points
    assert(out.pointlist != nullptr && "TETGENIO OUTPUT BUFFER EMPTY");
    v.resize(out.numberofpoints, std::vector<Real>(3));
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

bool MeshToTetgenio(tetgenio& in, const MatrixXr& V, const Eigen::MatrixXi& F) {

    std::vector<std::vector<Real>> v;
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

void Tetrahedralize(tetgenio& out, const MatrixXr& V, const Eigen::MatrixXi& F,
                    const std::string& flags) {
    tetgenio in;
    assert(MeshToTetgenio(in, V, F) && "FAILED TO CONVERT MESH TO TETGENIO");

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
