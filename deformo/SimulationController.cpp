#include "SimulationController.h"

SimulationController::SimulationController(const std::shared_ptr<Mesh>& mesh)
    : mesh_(mesh) {
    const auto boundary_conditions = GenerateDefaultBoundaryConditions(mesh_);
    simulation =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh_, boundary_conditions);
}

void SimulationController::StepForward() {
    simulation->Solve();
    recorded_displacements_.push_back(simulation->Integrate());
}

void SimulationController::Reset(const Ui::deformoClass& ui) {
    mesh_->Reset();
    const auto boundary_conditions = GenerateDefaultBoundaryConditions(mesh_);
    simulation =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh_, boundary_conditions);
    simulation->dt = dt_;
    simulation->current_time = 0.f;
    ResetPlaybackControls(ui);
}

void SimulationController::ResetPlaybackControls(const Ui::deformoClass& ui) {
    recorded_displacements_.clear();
    recorded_displacements_.resize(0);
    PlaybackPauseButtonPressed(ui);
}

void SimulationController::SetNodalMass(Real value) { nodal_mass_ = value; }

void SimulationController::SetPoissonsRatio(double value) {
    poissons_ratio_ = value;
}

void SimulationController::SetYoungsModulus(double value) {
    youngs_modulus_ = value;
}

void SimulationController::SetTimestepSize(double value) { dt_ = value; }

void SimulationController::SetRayleighLambda(double value) {
    rayleigh_lambda_ = value;
}

void SimulationController::SetRayleighMu(double value) { rayleigh_mu_ = value; }

void SimulationController::RunSimulationButtonPressed(
    const Ui::deformoClass& ui) {
    Reset(ui);
}

void SimulationController::PlaybackSkipStartButtonPressed(
    const Ui::deformoClass& ui) {
    ui.playback_controller->setValue(0);
}

void SimulationController::PlaybackSkipEndButtonPressed(
    const Ui::deformoClass& ui) {
    if (!recorded_displacements_.empty()) {
        ui.playback_controller->setValue(recorded_displacements_.size() - 1);
    }
}

void SimulationController::PlaybackPauseButtonPressed(
    const Ui::deformoClass& ui) {
    if (!recorded_displacements_.empty()) {
        ui.playback_controller->setMaximum(recorded_displacements_.size() - 1);
    }
    simulating = false;
}

void SimulationController::PlaybackPlayButtonPressed() { simulating = true; }

void SimulationController::PlaybackSliderChanged(int value) {
    mesh_->Update(recorded_displacements_.at(value));
}

void SimulationController::TabBarClicked(int index,
                                         const Ui::deformoClass& ui) {
    if (index == TabWindow::kDesigner) {
        if (!recorded_displacements_.empty()) {
            PlaybackPauseButtonPressed(ui);
            PlaybackSliderChanged(0);
        }
    }
}

BoundaryConditions SimulationController::GenerateDefaultBoundaryConditions(
    const std::shared_ptr<Mesh>& mesh) {
    const Eigen::Vector3f force(0.f, -100.f, 0.f);
    std::vector<unsigned int> indices;
    utils::FindMaxVertices(indices, mesh->positions);
    return AssignBoundaryConditionToFixedNodes(indices, force);
}
