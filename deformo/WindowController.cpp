#include "WindowController.h"
#include "Mesh.h"
#include "QTUtils.h"
#include "Utils.h"

WindowController::WindowController(Ui::deformoClass ui,
                                   const std::string& mesh_path) {
    const std::string suffix = ".ply";
    // TODO (@jparr721) - Make this a ui error
    assert(0 == mesh_path.compare(mesh_path.size() - suffix.size(),
                                  suffix.size(), suffix) &&
           "YOU CAN ONLY LOAD .PLY FILES");

    const auto mesh = std::make_shared<Mesh>(mesh_path, 0.f);
    const auto boundary_conditions = GenerateDefaultBoundaryConditions(mesh);
    simulation_ =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh, boundary_conditions);
}

void WindowController::SetSliceAxis(const QString& value) {
    slice_axis_ = utils::qt::QStringToString(value);
    emit OnSliceAxisChange(value);
}

void WindowController::SetSliceValue(float value) {
    slice_value_ = value;
    emit OnSliceValueChange(value);
}

void WindowController::SetNodalMass(float value) {
    nodal_mass_ = value;
    emit OnNodalMassChange(value);
}

void WindowController::SetPoissonsRatio(double value) {
    poissons_ratio_ = value;
    emit OnPoissonsRatioChange(value);
}

void WindowController::SetYoungsModulus(double value) {
    youngs_modulus_ = value;
    emit OnYoungsModulusChange(value);
}

void WindowController::SetTimestepSize(double value) {
    dt_ = value;
    emit OnTimestepSizeChange(value);
}

void WindowController::SetRayleighLambda(double value) {
    rayleigh_lambda_ = value;
    emit OnRayleighLambdaChange(value);
}

void WindowController::SetRayleighMu(double value) {
    rayleigh_mu_ = value;
    emit OnRayleighMuChange(value);
}

void WindowController::RunSimulationButtonPressed() {
    std::cout << "Run Simulation Pressed" << std::endl;
}

void WindowController::SetTetgenFlags(const QString& value) {
    tetgen_flags_ = utils::qt::QStringToString(value);
    emit OnTetgenFlagsChange(value);
}

void WindowController::RenderSimulationButtonPressed() {
    std::cout << "Render Button Pressed" << std::endl;
}

BoundaryConditions WindowController::GenerateDefaultBoundaryConditions(
    const std::shared_ptr<Mesh>& mesh) {
    const Eigen::Vector3f force(0.f, -9.8f, 0.f);

    std::vector<unsigned int> indices;
    utils::FindMaxVertices(indices, mesh->positions);
    return AssignBoundaryConditionToFixedNodes(indices, force);
}
