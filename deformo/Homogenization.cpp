#include "Homogenization.h"
#include "Utils.h"

Homogenization::Homogenization(std::shared_ptr<Rve> rve) {
    rve_ = rve;
    voxel_ = rve_->ToImplicitSurface();
    cell_len_x_ = rve_->width;
    cell_len_y_ = rve_->height;
    cell_len_z_ = rve_->depth;

    // For two-material composites, we sum the material parameters
    Tensor3r material_one_lambda = Where(voxel_, rve_->material_1.number);
    Tensor3r scalar_tensor_placeholder(material_one_lambda.Dimension(0),
                                       material_one_lambda.Dimension(1),
                                       material_one_lambda.Dimension(2));
    scalar_tensor_placeholder.SetConstant(rve_->material_1.Lambda());
    material_one_lambda.Instance() *= scalar_tensor_placeholder.Instance();
    Tensor3r material_two_lambda = Where(voxel_, rve_->material_2.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_2.Lambda());
    material_two_lambda.Instance() *= scalar_tensor_placeholder.Instance();

    lambda_ = Tensor3r(material_one_lambda.Instance() +
                       material_two_lambda.Instance());

    Tensor3r material_one_mu = Where(voxel_, rve_->material_1.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_1.Mu());
    material_one_mu.Instance() *= scalar_tensor_placeholder.Instance();
    Tensor3r material_two_mu = Where(voxel_, rve_->material_2.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_2.Mu());
    material_two_mu.Instance() *= scalar_tensor_placeholder.Instance();

    mu_ = Tensor3r(material_one_mu.Instance() + material_two_mu.Instance());
}

auto Homogenization::ComputeHexahedron(Real a, Real b, Real c)
    -> std::array<MatrixXr, 4> {
    // Constitutive matrix contribution for Mu
    Matrix6r C_mu = Matrix6r::Identity();
    C_mu(0, 0) = 2;
    C_mu(1, 1) = 2;
    C_mu(2, 2) = 2;

    // Constitutive matrix constribution for Lambda
    Matrix6r C_lambda;
    C_lambda.setZero();
    C_lambda(0, 0) = 1;
    C_lambda(1, 0) = 1;
    C_lambda(2, 0) = 1;
    C_lambda(0, 1) = 1;
    C_lambda(0, 2) = 1;
    C_lambda(1, 2) = 1;

    const auto xx = Vector3r(-std::sqrt(3. / 5.), 0, std::sqrt(3. / 5.));
    const auto yy = Vector3r(-std::sqrt(3. / 5.), 0, std::sqrt(3. / 5.));
    const auto zz = Vector3r(-std::sqrt(3. / 5.), 0, std::sqrt(3. / 5.));
    const auto ww = Vector3r(5. / 9., 8. / 9., 5. / 9.);

    MatrixXr ke_lambda = MatrixXr::Zero(24, 24);
    MatrixXr fe_lambda = MatrixXr::Zero(24, 6);

    MatrixXr ke_mu = MatrixXr::Zero(24, 24);
    MatrixXr fe_mu = MatrixXr::Zero(24, 6);

    for (int ii = 0; ii < xx.rows(); ++ii) {
        for (int jj = 0; jj < yy.rows(); ++jj) {
            for (int kk = 0; kk < zz.rows(); ++kk) {
                // Integration point
                const Real x = xx[ii];
                const Real y = yy[jj];
                const Real z = zz[kk];

                VectorXr qx;
                qx.resize(8);
                qx.row(0) << -((y - 1) * (z - 1)) / 8;
                qx.row(1) << ((y - 1) * (z - 1)) / 8;
                qx.row(2) << -((y + 1) * (z - 1)) / 8;
                qx.row(3) << ((y + 1) * (z - 1)) / 8;
                qx.row(4) << ((y - 1) * (z + 1)) / 8;
                qx.row(5) << -((y - 1) * (z + 1)) / 8;
                qx.row(6) << ((y + 1) * (z + 1)) / 8;
                qx.row(7) << -((y + 1) * (z + 1)) / 8;

                VectorXr qy;
                qy.resize(8);
                qy.row(0) << -((x - 1) * (z - 1)) / 8;
                qy.row(1) << ((x + 1) * (z - 1)) / 8;
                qy.row(2) << -((x + 1) * (z - 1)) / 8;
                qy.row(3) << ((x - 1) * (z - 1)) / 8;
                qy.row(4) << ((x - 1) * (z + 1)) / 8;
                qy.row(5) << -((x + 1) * (z + 1)) / 8;
                qy.row(6) << ((x + 1) * (z + 1)) / 8;
                qy.row(7) << -((x - 1) * (z + 1)) / 8;

                VectorXr qz;
                qz.resize(8);
                qz.row(0) << -((x - 1) * (y - 1)) / 8;
                qz.row(1) << ((x + 1) * (y - 1)) / 8;
                qz.row(2) << -((x + 1) * (y + 1)) / 8;
                qz.row(3) << ((x - 1) * (y + 1)) / 8;
                qz.row(4) << ((x - 1) * (y - 1)) / 8;
                qz.row(5) << -((x + 1) * (y - 1)) / 8;
                qz.row(6) << ((x + 1) * (y + 1)) / 8;
                qz.row(7) << -((x - 1) * (y + 1)) / 8;

                MatrixXr qq;
                qq.resize(3, 8);
                qq.row(0) = qx;
                qq.row(1) = qy;
                qq.row(2) = qz;

                MatrixXr dims;
                dims.resize(3, 8);
                dims.row(0) << -a, a, a, -a, -a, a, a, -a;
                dims.row(1) << -b, -b, b, b, -b, -b, b, b;
                dims.row(2) << -c, -c, -c, -c, c, c, c, c;
                dims.transposeInPlace();

                // Compute the jacobian matrix
                const MatrixXr J = qq * dims;
                const MatrixXr qxyz = J.fullPivLu().solve(qq);

                Tensor3r B_e = Tensor3r(8, 6, 3);
                B_e.SetConstant(0);
                const auto layers = B_e.Dimension(0);

                for (int layer = 0; layer < layers; ++layer) {
                    B_e(qxyz(0, layer), layer, 0, 0);

                    B_e(qxyz(1, layer), layer, 1, 1);

                    B_e(qxyz(2, layer), layer, 2, 2);

                    B_e(qxyz(1, layer), layer, 3, 0);
                    B_e(qxyz(0, layer), layer, 3, 1);

                    B_e(qxyz(2, layer), layer, 4, 1);
                    B_e(qxyz(1, layer), layer, 4, 2);

                    B_e(qxyz(2, layer), layer, 5, 0);
                    B_e(qxyz(0, layer), layer, 5, 2);
                }

                MatrixXr B = MatrixXr::Zero(6, 24);
                B << B_e.At(0), B_e.At(1), B_e.At(2), B_e.At(3), B_e.At(4),
                    B_e.At(5), B_e.At(6), B_e.At(7);
                const MatrixXr BT = B.transpose();

                const Real weight = J.determinant() * ww(ii) * ww(jj) * ww(kk);

                // Element stiffness coefficient matrices
                ke_lambda += weight * ((BT * C_lambda) * B);
                ke_mu += weight * ((BT * C_mu) * B);

                // Element load coefficient matrices
                fe_lambda += weight * (BT * C_lambda);
                fe_mu += weight * (BT * C_mu);
            }
        }
    }

    return std::array{ke_lambda, ke_mu, fe_lambda, fe_mu};
}
