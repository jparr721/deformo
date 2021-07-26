#include "Homogenization.h"

Homogenization::Homogenization(std::shared_ptr<Rve> rve) {
    rve_ = rve;
    voxel_ = rve_->ToImplicitSurface();
    cell_len_x_ = rve_->width;
    cell_len_y_ = rve_->height;
    cell_len_z_ = rve_->depth;

    // For two-material composites, we sum the material parameters
    Tensor3r material_one_lambda = Where(voxel_, rve_->material_1.number);
    material_one_lambda.Instance() *= rve_->material_1.Lambda();
    Tensor3r material_two_lambda = Where(voxel_, rve_->material_2.number);
    material_two_lambda.Instance() *= rve_->material_2.Lambda();

    lambda_ = Tensor3r(material_one_lambda.Instance() + material_two_lambda.Instance());

    Tensor3r material_one_mu = Where(voxel_, rve_->material_1.number);
    material_one_mu.Instance() *= rve_->material_1.Mu();
    Tensor3r material_two_mu = Where(voxel_, rve_->material_2.number);
    material_two_mu.Instance() *= rve_->material_2.Mu();

    mu_ = Tensor3r(material_one_mu.Instance() + material_two_mu.Instance());
}

auto Homogenization::ComputeHexahedron(Real a, Real b, Real c) -> Vector4r {
    // Constitutive matrix contribution for Mu
    Matrix6r C_mu = Matrix6r::Identity();
    C_mu(0, 0) = 2;
    C_mu(1, 1) = 2;
    C_mu(2, 2) = 2;

    // Constitutive matrix constribution for Lambda
    Matrix6r C_lambda;
    C_lambda(0, 0) = 1;
    C_lambda(1, 0) = 1;
    C_lambda(2, 0) = 1;
    C_lambda(0, 1) = 1;
    C_lambda(0, 2) = 1;
    C_lambda(1, 2) = 1;
}
