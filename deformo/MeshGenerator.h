#pragma once

#define TETLIBRARY

#include "Numerics.h"
#include <tetgen.h>

void Tetrahedralize(tetgenio& out, const MatrixXr& V, const Eigen::MatrixXi& F,
                    const std::string& flags);

bool MeshToTetgenio(tetgenio& in, const MatrixXr& V, const Eigen::MatrixXi& F);

bool TetgenioToMesh(MatrixXr& V, Eigen::MatrixXi& F, Eigen::MatrixXi& T,
                    const tetgenio& out);
