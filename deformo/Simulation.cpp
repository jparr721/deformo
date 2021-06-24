#include "Simulation.h"

Simulation::Simulation(
    float youngs_modulus, float poissons_ratio, float point_mass,
    const std::shared_ptr<Mesh>& mesh,
    const std::vector<BoundaryCondition>& boundary_conditions)
    : mesh_(mesh), youngs_modulus_(youngs_modulus),
      poissons_ratio_(poissons_ratio) {
    engine_ = std::make_unique<LinearTetrahedral>(
        youngs_modulus_, poissons_ratio_, mesh, boundary_conditions);
    integrator_ = std::make_unique<ExplicitCentralDifferenceMethod>(
        dt, point_mass, engine_->per_element_stiffness,
        engine_->global_displacement, engine_->boundary_forces);
}

void Simulation::Solve() {
    assert(engine_ != nullptr);
    assert(integrator_ != nullptr);
    current_time += dt;
    engine_->Solve(youngs_modulus_, poissons_ratio_, mesh_);
}

void Simulation::Integrate() const {
    integrator_->Solve(engine_->global_displacement, engine_->boundary_forces);
    mesh_->Update(engine_->ComputeRenderedDisplacements(mesh_->Size()));
}
