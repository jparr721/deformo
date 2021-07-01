#include "Simulation.h"

Simulation::Simulation(
    Real youngs_modulus, Real poissons_ratio, Real point_mass,
    const std::shared_ptr<Mesh>& mesh,
    const std::vector<BoundaryCondition>& boundary_conditions)
    : youngs_modulus_(youngs_modulus), poissons_ratio_(poissons_ratio),
      mesh_(mesh) {
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

VectorXr Simulation::Integrate() const {
    integrator_->Solve(engine_->global_displacement, engine_->boundary_forces);
    const VectorXr displacements =
        engine_->ComputeRenderedDisplacements(mesh_->Size());
    mesh_->Update(displacements);
    return displacements;
}
