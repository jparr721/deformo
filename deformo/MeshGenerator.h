#pragma once

#define TETLIBRARY

#include <Eigen/Dense>
#include <tetgen.h>

void Tetrahedralize(tetgenio& out, const Eigen::MatrixXf& V,
                    const Eigen::MatrixXi& F, const std::string& flags);

bool MeshToTetgenio(tetgenio& in, const Eigen::MatrixXf& V,
                    const Eigen::MatrixXi& F);

bool TetgenioToMesh(Eigen::MatrixXf& V, Eigen::MatrixXi& F, Eigen::MatrixXi& T,
                    const tetgenio& out);
