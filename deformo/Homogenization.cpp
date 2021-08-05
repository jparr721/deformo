#include "Homogenization.h"
#include "Utils.h"
#include <unsupported/Eigen/KroneckerProduct>

Homogenization::Homogenization(std::shared_ptr<Rve> rve) : rve_(rve) {
    // TODO(@jparr721) Remove (true) when done testing.
    voxel_ = rve_->ToImplicitSurface(true);
    cell_len_x_ = rve_->width;
    cell_len_y_ = rve_->height;
    cell_len_z_ = rve_->depth;

    // For two-material composites, we sum the material parameters
    Tensor3r material_one_lambda = voxel_.Where(rve_->material_1.number);
    Tensor3r scalar_tensor_placeholder(material_one_lambda.Dimensions());
    scalar_tensor_placeholder.SetConstant(rve_->material_1.lambda);
    material_one_lambda.Instance() *= scalar_tensor_placeholder.Instance();

    Tensor3r material_two_lambda = voxel_.Where(rve_->material_2.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_2.lambda);
    material_two_lambda.Instance() *= scalar_tensor_placeholder.Instance();

    lambda_ = Tensor3r(material_one_lambda.Instance() +
                       material_two_lambda.Instance());

    Tensor3r material_one_mu = voxel_.Where(rve_->material_1.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_1.G);
    material_one_mu.Instance() *= scalar_tensor_placeholder.Instance();
    Tensor3r material_two_mu = voxel_.Where(rve_->material_2.number);
    scalar_tensor_placeholder.SetConstant(rve_->material_2.G);
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
                const MatrixXr qxyz = J.ldlt().solve(qq);

                Tensor3r B_e = Tensor3r(6, 3, 8);
                B_e.SetConstant(0);
                const auto layers = B_e.Dimension(2);

                for (int layer = 0; layer < layers; ++layer) {
                    B_e(0, 0, layer) = qxyz(0, layer);
                    B_e(1, 1, layer) = qxyz(1, layer);
                    B_e(2, 2, layer) = qxyz(2, layer);
                    B_e(3, 0, layer) = qxyz(1, layer);
                    B_e(3, 1, layer) = qxyz(0, layer);
                    B_e(4, 1, layer) = qxyz(2, layer);
                    B_e(4, 2, layer) = qxyz(1, layer);
                    B_e(5, 0, layer) = qxyz(2, layer);
                    B_e(5, 2, layer) = qxyz(0, layer);
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

auto Homogenization::ComputeElementDegreesOfFreedom(unsigned int n_elements)
    -> MatrixXi {
    assertion.Assert(voxel_.Dimensions().size() == 3, __FUNCTION__, __FILE__,
                     __LINE__, "Voxel is improperly shaped");
    const unsigned int n_el_x = voxel_.Dimension(0);
    const unsigned int n_el_y = voxel_.Dimension(1);
    const unsigned int n_el_z = voxel_.Dimension(2);

    const unsigned int number_of_nodes =
        (1 + n_el_x) * (1 + n_el_y) * (1 + n_el_z);

    // Set up to apply the periodic boundary conditions for periodic volumes.
    // Here, we set up the node numbers and indexing degrees of freedom for
    // 3-D Homogenization.
    VectorXi _nn = VectorXi::LinSpaced(number_of_nodes, 1, number_of_nodes);

    assertion.Assert(_nn.size() == number_of_nodes, __FUNCTION__, __FILE__,
                     __LINE__, "Node numbers improperly formatted!",
                     number_of_nodes);
    const Tensor3i node_numbers =
        Tensor3i::Expand(_nn, 1 + n_el_x, 1 + n_el_y, 1 + n_el_z);

    const unsigned int node_numbers_x = node_numbers.Dimension(0) - 1;
    const unsigned int node_numbers_y = node_numbers.Dimension(1) - 1;
    const unsigned int node_numbers_z = node_numbers.Dimension(2) - 1;

    Tensor3i _dof(node_numbers_x, node_numbers_y, node_numbers_z);
    for (auto x = 0u; x < node_numbers_x; ++x) {
        for (auto y = 0u; y < node_numbers_y; ++y) {
            for (auto z = 0u; z < node_numbers_z; ++z) {
                _dof(x, y, z) = node_numbers.At(x, y, z);
            }
        }
    }

    Tensor3i three(_dof.Dimensions());
    three.SetConstant(3);
    _dof.Instance() *= three.Instance();
    Tensor3i one(_dof.Dimensions());
    one.SetConstant(1);
    _dof.Instance() += one.Instance();

    const VectorXi degrees_of_freedom =
        _dof.Matrix(_dof.Dimensions().prod(), 1);

    Vector6<int> _mid;
    _mid << 3, 4, 5, 0, 1, 2;
    _mid += Vector6<int>::Ones() * 3 * n_el_x;

    Vector12<int> _add_x;
    _add_x << 0, 1, 2, _mid, -3, -2, -1;

    const Vector12<int> _add_xy =
        (_add_x.array() + (3 * (1 + n_el_y) * (1 + n_el_x))).matrix();

    MatrixXi _add_combined(1, 24);
    _add_combined << _add_x.transpose(), _add_xy.transpose();

    const MatrixXi _edof_lhs = degrees_of_freedom.replicate(1, 24);
    const MatrixXi _edof_rhs = _add_combined.replicate(n_elements, 1);

    return _edof_lhs + _edof_rhs;
}

auto Homogenization::ComputeUniqueNodes(unsigned int n_elements) -> Tensor3i {
    assertion.Assert(voxel_.Dimensions().size() == 3, __FUNCTION__, __FILE__,
                     __LINE__, "Voxel is improperly shaped");
    const unsigned int n_el_x = voxel_.Dimension(0);
    const unsigned int n_el_y = voxel_.Dimension(1);
    const unsigned int n_el_z = voxel_.Dimension(2);

    const VectorXi _uniq_el = VectorXi::LinSpaced(n_elements, 1, n_elements);

    Tensor3i _uniq_t_1 = Tensor3i::Expand(_uniq_el, n_el_x, n_el_y, n_el_z);

    Tensor3i _index_tensor((_uniq_t_1.Dimensions().array() + 1).matrix());

    // Extend with a mirror of the back border
    std::vector<VectorXi> back_borders;
    constexpr int row = 0;
    for (auto layer_idx = 0u; layer_idx < n_el_z; ++layer_idx) {
        back_borders.emplace_back(_uniq_t_1.At(layer_idx, row));
    }

    Tensor3i _uniq_t_2 =
        _uniq_t_1.Append(back_borders, Tensor3i::InsertOpIndex::kEnd,
                         Tensor3i::OpOrientation::kRow);

    // Extend with a mirror of the left border
    std::vector<VectorXi> left_borders;
    constexpr int col = 0;
    for (auto layer_idx = 0u; layer_idx < n_el_z; ++layer_idx) {
        left_borders.emplace_back(_uniq_t_2.Col(layer_idx, col));
    }

    Tensor3i _uniq_t_3 =
        _uniq_t_2.Append(left_borders, Tensor3i::InsertOpIndex::kEnd,
                         Tensor3i::OpOrientation::kCol);

    // Finally, extend with a mirror of the top border
    const MatrixXi first_layer = _uniq_t_3.At(0);
    return _uniq_t_3.Append(first_layer, Tensor3i::InsertOpIndex::kEnd);
}

auto Homogenization::ComputeUniqueDegreesOfFreedom(
    const MatrixXi& element_degrees_of_freedom, const Tensor3i& unique_nodes)
    -> MatrixXi {
    const unsigned int n_nodes =
        (voxel_.Dimensions().array() + 1).matrix().prod();

    VectorXi _dof = VectorXi::Ones(3 * n_nodes);
    VectorXi _uniq_vec = unique_nodes.Vector(unique_nodes.Dimensions().prod());

    for (int i = 0; i < _dof.rows(); i += 3) {
        const int idx = i / 3;
        _dof(i) = 3 * _uniq_vec(idx) - 2;
    }

    for (int i = 1; i < _dof.rows(); i += 3) {
        const int idx = i / 3;
        _dof(i) = 3 * _uniq_vec(idx) - 1;
    }

    for (int i = 2; i < _dof.rows(); i += 3) {
        const int idx = i / 3;
        _dof(i) = 3 * _uniq_vec(idx);
    }

    const MatrixXi indices = (element_degrees_of_freedom.array() - 1).matrix();

    return linear_algebra::IndexVectorByMatrix(_dof, indices);
}

auto Homogenization::AssembleStiffnessMatrix(
    const unsigned int n_degrees_of_freedom,
    const MatrixXi& unique_degrees_of_freedom, const MatrixXr& ke_lambda,
    const MatrixXr& ke_mu) -> SparseMatrixXr {
    const MatrixXi idx_i_kron =
        Eigen::kroneckerProduct(unique_degrees_of_freedom,
                                MatrixXi::Ones(24, 1))
            .adjoint();
    const VectorXi idx_i =
        ((linear_algebra::MatrixToVector(idx_i_kron)).array() - 1).matrix();

    const MatrixXi idx_j_kron =
        Eigen::kroneckerProduct(unique_degrees_of_freedom,
                                MatrixXi::Ones(1, 24))
            .adjoint();
    const VectorXi idx_j =
        ((linear_algebra::MatrixToVector(idx_j_kron)).array() - 1).matrix();

    const MatrixXr sK =
        (linear_algebra::MatrixToVector(ke_lambda) *
         lambda_.Vector(lambda_.Dimensions().prod()).transpose()) +
        (linear_algebra::MatrixToVector(ke_mu) *
         mu_.Vector(mu_.Dimensions().prod()).transpose());

    const VectorXr stiffness_entries = linear_algebra::MatrixToVector(sK);

    const std::vector<Eigen::Triplet<Real>> K_entries =
        linear_algebra::ToTriplets(idx_i, idx_j, stiffness_entries);

    SparseMatrixXr K(n_degrees_of_freedom, n_degrees_of_freedom);
    K.setFromTriplets(K_entries.begin(), K_entries.end());

    SparseMatrixXr KT = K.adjoint();

    K += KT;

    return K * 1/2;
}

auto Homogenization::AssembleLoadMatrix(
    const MatrixXi& unique_degrees_of_freedom) -> MatrixXr {
    return MatrixXr();
}
