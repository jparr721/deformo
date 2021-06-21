#define _USE_MATH_DEFINES

#include "LinearTetrahedral.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <utility>

#include "Rayleigh.h"
#include "Utils.h"

LinearTetrahedral::LinearTetrahedral(
    const float modulus_of_elasticity, const float poissons_ratio,
    const float point_mass, std::shared_ptr<Mesh> mesh,
    std::vector<BoundaryCondition> boundary_conditions)
    : modulus_of_elasticity(modulus_of_elasticity),
      poissons_ratio(poissons_ratio), mesh(std::move(mesh)),
      boundary_conditions(std::move(boundary_conditions)),
      integrator_size(boundary_conditions.size()) {
    AssembleElementStiffness();
    AssembleGlobalStiffness();
    AssembleMassMatrix(point_mass);
    AssembleBoundaryForces();
    InitializeIntegrator();
}

void LinearTetrahedral::Update() {
    current_time += dt;
    integrator->Solve(global_displacement, boundary_forces);
    mesh->Update(ComputeRenderedDisplacements());
}

void LinearTetrahedral::ComputeElementStiffness(
    Eigen::Matrix12f& element_stiffness, const Eigen::Vector3f& shape_one,
    const Eigen::Vector3f& shape_two, const Eigen::Vector3f& shape_three,
    const Eigen::Vector3f& shape_four) {
    const Eigen::MatrixXf B = AssembleStrainRelationshipMatrix(
        shape_one, shape_two, shape_three, shape_four);

    const Eigen::Matrix66f D =
        AssembleStressStrainMatrix(poissons_ratio, modulus_of_elasticity);

    const float V =
        ComputeElementVolume(shape_one, shape_two, shape_three, shape_four);
    element_stiffness = V * B.transpose() * D * B;
}

void LinearTetrahedral::AssembleElementStiffness() {
    for (int i = 0; i < mesh->SimNodesSize(); i += kFaceStride) {
        std::vector<int> stiffness_coordinates;
        // Get the index face value
        int index = mesh->GetPositionAtFaceIndex(i);
        stiffness_coordinates.push_back(index / 3);
        Eigen::VectorXf shape_one;
        utils::SliceEigenVector(shape_one, mesh->positions, index, index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 1);
        stiffness_coordinates.push_back(index / 3);
        Eigen::VectorXf shape_two;
        utils::SliceEigenVector(shape_two, mesh->positions, index, index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 2);
        stiffness_coordinates.push_back(index / 3);
        Eigen::VectorXf shape_three;
        utils::SliceEigenVector(shape_three, mesh->positions, index, index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 3);
        stiffness_coordinates.push_back(index / 3);
        Eigen::VectorXf shape_four;
        utils::SliceEigenVector(shape_four, mesh->positions, index, index + 2);

        Eigen::Matrix12f element_stiffness_matrix;
        ComputeElementStiffness(element_stiffness_matrix, shape_one, shape_two,
                                shape_three, shape_four);

        const ElementStiffness element_stiffness = {element_stiffness_matrix,
                                                    stiffness_coordinates};

        element_stiffnesses.emplace_back(element_stiffness);
    }
}

void LinearTetrahedral::AssembleElementStresses(const Eigen::VectorXf& u,
                                                const Eigen::MatrixXf& B) {
    const Eigen::Matrix66f D =
        AssembleStressStrainMatrix(poissons_ratio, modulus_of_elasticity);

    const Eigen::Vector6f sigma = D * B * u;

    sigmas.emplace_back(sigma);
}

void LinearTetrahedral::AssembleBoundaryForces() {
    assert(!boundary_conditions.empty() && "NO CONDITIONS TO SOLVE FOR");

    boundary_forces.resize(integrator_size * 3);
    boundary_forces.setZero();

    boundary_force_indices.resize(integrator_size * 3);

    int segment = 0;
    for (const auto& [node, force] : boundary_conditions) {
        // Get the node number so we can begin indexing
        const unsigned int node_number = node;

        // The row is the same as the index segment
        boundary_forces.segment(segment, 3) << force;

        boundary_force_indices.segment(segment, 3) << node_number,
            node_number + 1, node_number + 2;
        segment += 3;
    }

    utils::SliceByIndices(per_element_stiffness, global_stiffness,
                          boundary_force_indices, boundary_force_indices);
}

void LinearTetrahedral::AssembleElementPlaneStresses() {
    for (const auto& sigma : sigmas) {
        const float s1 = sigma.sum();
        const float s2 =
            (sigma(0) * sigma(1) + sigma(0) * sigma(2) + sigma(1) * sigma(2)) -
            (sigma(3) * sigma(3) - sigma(4) * sigma(4) - sigma(5) * sigma(5));

        Eigen::Matrix3f ms3;
        ms3.row(0) << sigma(0), sigma(3), sigma(5);
        ms3.row(1) << sigma(3), sigma(1), sigma(4);
        ms3.row(2) << sigma(5), sigma(4), sigma(2);

        const float s3 = ms3.determinant();

        const Eigen::Vector3f plane_stress(s1, s2, s3);
        element_stresses.emplace_back(plane_stress);
    }
}

void LinearTetrahedral::AssembleGlobalStiffness() {
    // Allocate space in global_stiffness for 3nx3n elements
    const unsigned int size = mesh->Size();

    global_stiffness = Eigen::MatrixXf::Zero(size, size);

    for (const auto& element_stiffness : element_stiffnesses) {
        const Eigen::Matrix12f k = element_stiffness.stiffness_matrix;
        const auto i = element_stiffness.indices.at(0);
        const auto j = element_stiffness.indices.at(1);
        const auto m = element_stiffness.indices.at(2);
        const auto n = element_stiffness.indices.at(3);

        global_stiffness(3 * i, 3 * i) += k(0, 0);
        global_stiffness(3 * i, 3 * i + 1) += k(0, 1);
        global_stiffness(3 * i, 3 * i + 2) += k(0, 2);
        global_stiffness(3 * i, 3 * j) += k(0, 3);
        global_stiffness(3 * i, 3 * j + 1) += k(0, 4);
        global_stiffness(3 * i, 3 * j + 2) += k(0, 5);
        global_stiffness(3 * i, 3 * m) += k(0, 6);
        global_stiffness(3 * i, 3 * m + 1) += k(0, 7);
        global_stiffness(3 * i, 3 * m + 2) += k(0, 8);
        global_stiffness(3 * i, 3 * n) += k(0, 9);
        global_stiffness(3 * i, 3 * n + 1) += k(0, 10);
        global_stiffness(3 * i, 3 * n + 2) += k(0, 11);

        global_stiffness(3 * i + 1, 3 * i) += k(1, 0);
        global_stiffness(3 * i + 1, 3 * i + 1) += k(1, 1);
        global_stiffness(3 * i + 1, 3 * i + 2) += k(1, 2);
        global_stiffness(3 * i + 1, 3 * j) += k(1, 3);
        global_stiffness(3 * i + 1, 3 * j + 1) += k(1, 4);
        global_stiffness(3 * i + 1, 3 * j + 2) += k(1, 5);
        global_stiffness(3 * i + 1, 3 * m) += k(1, 6);
        global_stiffness(3 * i + 1, 3 * m + 1) += k(1, 7);
        global_stiffness(3 * i + 1, 3 * m + 2) += k(1, 8);
        global_stiffness(3 * i + 1, 3 * n) += k(1, 9);
        global_stiffness(3 * i + 1, 3 * n + 1) += k(1, 10);
        global_stiffness(3 * i + 1, 3 * n + 2) += k(1, 11);

        global_stiffness(3 * i + 2, 3 * i) += k(2, 0);
        global_stiffness(3 * i + 2, 3 * i + 1) += k(2, 1);
        global_stiffness(3 * i + 2, 3 * i + 2) += k(2, 2);
        global_stiffness(3 * i + 2, 3 * j) += k(2, 3);
        global_stiffness(3 * i + 2, 3 * j + 1) += k(2, 4);
        global_stiffness(3 * i + 2, 3 * j + 2) += k(2, 5);
        global_stiffness(3 * i + 2, 3 * m) += k(2, 6);
        global_stiffness(3 * i + 2, 3 * m + 1) += k(2, 7);
        global_stiffness(3 * i + 2, 3 * m + 2) += k(2, 8);
        global_stiffness(3 * i + 2, 3 * n) += k(2, 9);
        global_stiffness(3 * i + 2, 3 * n + 1) += k(2, 10);
        global_stiffness(3 * i + 2, 3 * n + 2) += k(2, 11);

        // j
        global_stiffness(3 * j, 3 * i) += k(3, 0);
        global_stiffness(3 * j, 3 * i + 1) += k(3, 1);
        global_stiffness(3 * j, 3 * i + 2) += k(3, 2);
        global_stiffness(3 * j, 3 * j) += k(3, 3);
        global_stiffness(3 * j, 3 * j + 1) += k(3, 4);
        global_stiffness(3 * j, 3 * j + 2) += k(3, 5);
        global_stiffness(3 * j, 3 * m) += k(3, 6);
        global_stiffness(3 * j, 3 * m + 1) += k(3, 7);
        global_stiffness(3 * j, 3 * m + 2) += k(3, 8);
        global_stiffness(3 * j, 3 * n) += k(3, 9);
        global_stiffness(3 * j, 3 * n + 1) += k(3, 10);
        global_stiffness(3 * j, 3 * n + 2) += k(3, 11);

        global_stiffness(3 * j + 1, 3 * i) += k(4, 0);
        global_stiffness(3 * j + 1, 3 * i + 1) += k(4, 1);
        global_stiffness(3 * j + 1, 3 * i + 2) += k(4, 2);
        global_stiffness(3 * j + 1, 3 * j) += k(4, 3);
        global_stiffness(3 * j + 1, 3 * j + 1) += k(4, 4);
        global_stiffness(3 * j + 1, 3 * j + 2) += k(4, 5);
        global_stiffness(3 * j + 1, 3 * m) += k(4, 6);
        global_stiffness(3 * j + 1, 3 * m + 1) += k(4, 7);
        global_stiffness(3 * j + 1, 3 * m + 2) += k(4, 8);
        global_stiffness(3 * j + 1, 3 * n) += k(4, 9);
        global_stiffness(3 * j + 1, 3 * n + 1) += k(4, 10);
        global_stiffness(3 * j + 1, 3 * n + 2) += k(4, 11);

        global_stiffness(3 * j + 2, 3 * i) += k(5, 0);
        global_stiffness(3 * j + 2, 3 * i + 1) += k(5, 1);
        global_stiffness(3 * j + 2, 3 * i + 2) += k(5, 2);
        global_stiffness(3 * j + 2, 3 * j) += k(5, 3);
        global_stiffness(3 * j + 2, 3 * j + 1) += k(5, 4);
        global_stiffness(3 * j + 2, 3 * j + 2) += k(5, 5);
        global_stiffness(3 * j + 2, 3 * m) += k(5, 6);
        global_stiffness(3 * j + 2, 3 * m + 1) += k(5, 7);
        global_stiffness(3 * j + 2, 3 * m + 2) += k(5, 8);
        global_stiffness(3 * j + 2, 3 * n) += k(5, 9);
        global_stiffness(3 * j + 2, 3 * n + 1) += k(5, 10);
        global_stiffness(3 * j + 2, 3 * n + 2) += k(5, 11);

        // m
        global_stiffness(3 * m, 3 * i) += k(6, 0);
        global_stiffness(3 * m, 3 * i + 1) += k(6, 1);
        global_stiffness(3 * m, 3 * i + 2) += k(6, 2);
        global_stiffness(3 * m, 3 * j) += k(6, 3);
        global_stiffness(3 * m, 3 * j + 1) += k(6, 4);
        global_stiffness(3 * m, 3 * j + 2) += k(6, 5);
        global_stiffness(3 * m, 3 * m) += k(6, 6);
        global_stiffness(3 * m, 3 * m + 1) += k(6, 7);
        global_stiffness(3 * m, 3 * m + 2) += k(6, 8);
        global_stiffness(3 * m, 3 * n) += k(6, 9);
        global_stiffness(3 * m, 3 * n + 1) += k(6, 10);
        global_stiffness(3 * m, 3 * n + 2) += k(6, 11);

        global_stiffness(3 * m + 1, 3 * i) += k(7, 0);
        global_stiffness(3 * m + 1, 3 * i + 1) += k(7, 1);
        global_stiffness(3 * m + 1, 3 * i + 2) += k(7, 2);
        global_stiffness(3 * m + 1, 3 * j) += k(7, 3);
        global_stiffness(3 * m + 1, 3 * j + 1) += k(7, 4);
        global_stiffness(3 * m + 1, 3 * j + 2) += k(7, 5);
        global_stiffness(3 * m + 1, 3 * m) += k(7, 6);
        global_stiffness(3 * m + 1, 3 * m + 1) += k(7, 7);
        global_stiffness(3 * m + 1, 3 * m + 2) += k(7, 8);
        global_stiffness(3 * m + 1, 3 * n) += k(7, 9);
        global_stiffness(3 * m + 1, 3 * n + 1) += k(7, 10);
        global_stiffness(3 * m + 1, 3 * n + 2) += k(7, 11);

        global_stiffness(3 * m + 2, 3 * i) += k(8, 0);
        global_stiffness(3 * m + 2, 3 * i + 1) += k(8, 1);
        global_stiffness(3 * m + 2, 3 * i + 2) += k(8, 2);
        global_stiffness(3 * m + 2, 3 * j) += k(8, 3);
        global_stiffness(3 * m + 2, 3 * j + 1) += k(8, 4);
        global_stiffness(3 * m + 2, 3 * j + 2) += k(8, 5);
        global_stiffness(3 * m + 2, 3 * m) += k(8, 6);
        global_stiffness(3 * m + 2, 3 * m + 1) += k(8, 7);
        global_stiffness(3 * m + 2, 3 * m + 2) += k(8, 8);
        global_stiffness(3 * m + 2, 3 * n) += k(8, 9);
        global_stiffness(3 * m + 2, 3 * n + 1) += k(8, 10);
        global_stiffness(3 * m + 2, 3 * n + 2) += k(8, 11);

        // n
        global_stiffness(3 * n, 3 * i) += k(9, 0);
        global_stiffness(3 * n, 3 * i + 1) += k(9, 1);
        global_stiffness(3 * n, 3 * i + 2) += k(9, 2);
        global_stiffness(3 * n, 3 * j) += k(9, 3);
        global_stiffness(3 * n, 3 * j + 1) += k(9, 4);
        global_stiffness(3 * n, 3 * j + 2) += k(9, 5);
        global_stiffness(3 * n, 3 * m) += k(9, 6);
        global_stiffness(3 * n, 3 * m + 1) += k(9, 7);
        global_stiffness(3 * n, 3 * m + 2) += k(9, 8);
        global_stiffness(3 * n, 3 * n) += k(9, 9);
        global_stiffness(3 * n, 3 * n + 1) += k(9, 10);
        global_stiffness(3 * n, 3 * n + 2) += k(9, 11);

        global_stiffness(3 * n + 1, 3 * i) += k(10, 0);
        global_stiffness(3 * n + 1, 3 * i + 1) += k(10, 1);
        global_stiffness(3 * n + 1, 3 * i + 2) += k(10, 2);
        global_stiffness(3 * n + 1, 3 * j) += k(10, 3);
        global_stiffness(3 * n + 1, 3 * j + 1) += k(10, 4);
        global_stiffness(3 * n + 1, 3 * j + 2) += k(10, 5);
        global_stiffness(3 * n + 1, 3 * m) += k(10, 6);
        global_stiffness(3 * n + 1, 3 * m + 1) += k(10, 7);
        global_stiffness(3 * n + 1, 3 * m + 2) += k(10, 8);
        global_stiffness(3 * n + 1, 3 * n) += k(10, 9);
        global_stiffness(3 * n + 1, 3 * n + 1) += k(10, 10);
        global_stiffness(3 * n + 1, 3 * n + 2) += k(10, 11);

        global_stiffness(3 * n + 2, 3 * i) += k(11, 0);
        global_stiffness(3 * n + 2, 3 * i + 1) += k(11, 1);
        global_stiffness(3 * n + 2, 3 * i + 2) += k(11, 2);
        global_stiffness(3 * n + 2, 3 * j) += k(11, 3);
        global_stiffness(3 * n + 2, 3 * j + 1) += k(11, 4);
        global_stiffness(3 * n + 2, 3 * j + 2) += k(11, 5);
        global_stiffness(3 * n + 2, 3 * m) += k(11, 6);
        global_stiffness(3 * n + 2, 3 * m + 1) += k(11, 7);
        global_stiffness(3 * n + 2, 3 * m + 2) += k(11, 8);
        global_stiffness(3 * n + 2, 3 * n) += k(11, 9);
        global_stiffness(3 * n + 2, 3 * n + 1) += k(11, 10);
        global_stiffness(3 * n + 2, 3 * n + 2) += k(11, 11);
    }
}

Eigen::Matrix66f
LinearTetrahedral::AssembleStressStrainMatrix(float poissons_ratio,
                                              float modulus_of_elasticity) {
    Eigen::Matrix66f D;
    D.row(0) << 1 - poissons_ratio, poissons_ratio, poissons_ratio, 0, 0, 0;
    D.row(1) << poissons_ratio, 1 - poissons_ratio, poissons_ratio, 0, 0, 0;
    D.row(2) << poissons_ratio, poissons_ratio, 1 - poissons_ratio, 0, 0, 0;
    D.row(3) << 0, 0, 0, (1 - 2 * poissons_ratio) / 2, 0, 0;
    D.row(4) << 0, 0, 0, 0, (1 - 2 * poissons_ratio) / 2, 0;
    D.row(5) << 0, 0, 0, 0, 0, (1 - 2 * poissons_ratio) / 2;
    D *= modulus_of_elasticity /
         ((1 + poissons_ratio) * (1 - 2 * poissons_ratio));
    return D;
}

Eigen::MatrixXf LinearTetrahedral::AssembleStrainRelationshipMatrix(
    const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
    const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four) {
    Eigen::MatrixXf strain_relationship;
    const float V =
        ComputeElementVolume(shape_one, shape_two, shape_three, shape_four);
    const auto create_beta_submatrix = [](float beta, float gamma,
                                          float delta) -> BetaSubMatrixXf {
        BetaSubMatrixXf B;
        B.row(0) << beta, 0, 0;
        B.row(1) << 0, gamma, 0;
        B.row(2) << 0, 0, delta;
        B.row(3) << gamma, beta, 0;
        B.row(4) << 0, delta, gamma;
        B.row(5) << delta, 0, beta;
        return B;
    };

    const float x1 = shape_one.x();
    const float y1 = shape_one.y();
    const float z1 = shape_one.z();

    const float x2 = shape_two.x();
    const float y2 = shape_two.y();
    const float z2 = shape_two.z();

    const float x3 = shape_three.x();
    const float y3 = shape_three.y();
    const float z3 = shape_three.z();

    const float x4 = shape_four.x();
    const float y4 = shape_four.y();
    const float z4 = shape_four.z();

    const float beta_1 =
        -1 * ConstructShapeFunctionParameter(y2, z2, y3, z3, y4, z4);
    const float beta_2 =
        ConstructShapeFunctionParameter(y1, z1, y3, z3, y4, z4);
    const float beta_3 =
        -1 * ConstructShapeFunctionParameter(y1, z1, y2, z2, y4, z4);
    const float beta_4 =
        ConstructShapeFunctionParameter(y1, z1, y2, z2, y3, z3);

    const float gamma_1 =
        ConstructShapeFunctionParameter(x2, z2, x3, z3, x4, z4);
    const float gamma_2 =
        -1 * ConstructShapeFunctionParameter(x1, z1, x3, z3, x4, z4);
    const float gamma_3 =
        ConstructShapeFunctionParameter(x1, z1, x2, z2, x4, z4);
    const float gamma_4 =
        -1 * ConstructShapeFunctionParameter(x1, z1, x2, z2, x3, z3);

    const float delta_1 =
        -1 * ConstructShapeFunctionParameter(x2, y2, x3, y3, x4, y4);
    const float delta_2 =
        ConstructShapeFunctionParameter(x1, y1, x3, y3, x4, y4);
    const float delta_3 =
        -1 * ConstructShapeFunctionParameter(x1, y1, x2, y2, x4, y4);
    const float delta_4 =
        ConstructShapeFunctionParameter(x1, y1, x2, y2, x3, y3);

    const BetaSubMatrixXf B1 = create_beta_submatrix(beta_1, gamma_1, delta_1);
    const BetaSubMatrixXf B2 = create_beta_submatrix(beta_2, gamma_2, delta_2);
    const BetaSubMatrixXf B3 = create_beta_submatrix(beta_3, gamma_3, delta_3);
    const BetaSubMatrixXf B4 = create_beta_submatrix(beta_4, gamma_4, delta_4);

    // Matrix is 6 x 12
    strain_relationship.resize(B1.rows(), B1.cols() * 4);
    strain_relationship << B1, B2, B3, B4;
    if (V != 0.f) {
        strain_relationship /= (6 * V);
    } else {
        strain_relationship /= (6);
    }
    return strain_relationship;
}

float LinearTetrahedral::ConstructShapeFunctionParameter(float p1, float p2,
                                                         float p3, float p4,
                                                         float p5, float p6) {
    Eigen::Matrix3f parameter;
    parameter.row(0) << 1, p1, p2;
    parameter.row(1) << 1, p3, p4;
    parameter.row(2) << 1, p5, p6;
    return parameter.determinant();
}

Eigen::VectorXf LinearTetrahedral::ComputeRenderedDisplacements() {
    assert(!boundary_conditions.empty() && "NO BOUNDARY CONDITIONS");
    Eigen::VectorXf output = Eigen::VectorXf::Zero(mesh->Size());

    int i = 0;
    for (const auto& [node, _] : boundary_conditions) {
        output.segment(node, 3) << global_displacement(i),
            global_displacement(i + 1), global_displacement(i + 2);
        i += 3;
    }

    return output;
}

void LinearTetrahedral::AssembleMassMatrix(const float point_mass) {
    mass.resize(integrator_size * 3, integrator_size * 3);
    mass.setIdentity();
}

void LinearTetrahedral::Solve() {
    const Eigen::VectorXf solved_displacement = ComputeRenderedDisplacements();
    for (int i = 0; i < mesh->SimNodesSize(); i += kFaceStride) {
        int index = mesh->GetPositionAtFaceIndex(i);
        Eigen::VectorXf shape_one;
        utils::SliceEigenVector(shape_one, mesh->positions, index, index + 2);
        Eigen::VectorXf displacement_one;
        utils::SliceEigenVector(displacement_one, solved_displacement, index,
                                index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 1);
        Eigen::VectorXf shape_two;
        utils::SliceEigenVector(shape_two, mesh->positions, index, index + 2);
        Eigen::VectorXf displacement_two;
        utils::SliceEigenVector(displacement_two, solved_displacement, index,
                                index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 2);
        Eigen::VectorXf shape_three;
        utils::SliceEigenVector(shape_three, mesh->positions, index, index + 2);
        Eigen::VectorXf displacement_three;
        utils::SliceEigenVector(displacement_three, solved_displacement, index,
                                index + 2);

        index = mesh->GetPositionAtFaceIndex(i + 3);
        Eigen::VectorXf shape_four;
        utils::SliceEigenVector(shape_four, mesh->positions, index, index + 2);
        Eigen::VectorXf displacement_four;
        utils::SliceEigenVector(displacement_four, solved_displacement, index,
                                index + 2);

        const Eigen::MatrixXf B = AssembleStrainRelationshipMatrix(
            shape_one, shape_two, shape_three, shape_four);

        Eigen::VectorXf u(12);
        u << displacement_one, displacement_two, displacement_three,
            displacement_four;
        AssembleElementStresses(u, B);
    }

    AssembleElementPlaneStresses();
    Update();
}

float LinearTetrahedral::ComputeElementVolume(
    const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
    const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four) {
    const float x1 = shape_one.x();
    const float y1 = shape_one.y();
    const float z1 = shape_one.z();

    const float x2 = shape_two.x();
    const float y2 = shape_two.y();
    const float z2 = shape_two.z();

    const float x3 = shape_three.x();
    const float y3 = shape_three.y();
    const float z3 = shape_three.z();

    const float x4 = shape_four.x();
    const float y4 = shape_four.y();
    const float z4 = shape_four.z();

    Eigen::Matrix4f V;
    V.row(0) << 1, x1, y1, z1;
    V.row(1) << 1, x2, y2, z2;
    V.row(2) << 1, x3, y3, z3;
    V.row(3) << 1, x4, y4, z4;

    return V.determinant() / 6;
}

void LinearTetrahedral::InitializeIntegrator() {
    global_displacement.resize(integrator_size * 3);
    global_displacement.setZero();

    Eigen::MatrixXf damping;
    ComputeRayleighDamping(damping, per_element_stiffness, mass, 0.5, 0.5, 1e1);

    integrator = std::make_unique<ExplicitCentralDifferenceMethod>(
        dt, global_displacement, per_element_stiffness, mass, boundary_forces,
        damping);
}
