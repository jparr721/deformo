#pragma once

#include "Numerics.h"

void ComputeRayleighDamping(MatrixXr& out, const MatrixXr& stiffness,
                            const MatrixXr& mass, Real mu, Real lambda,
                            Real modifier = 1.f);
