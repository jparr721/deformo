#pragma once

#include "BoundaryCondition.h"
#include "ExplicitCentralDifference.h"
#include "LinearTetrahedral.h"
#include "Mesh.h"

#include <memory>
#include <vector>

class Simulation {
  public:
    float dt = 0.1f;
    float current_time = 0.f;

    Simulation() = default;
    Simulation(float youngs_modulus, float poissons_ratio, float point_mass,
               const std::shared_ptr<Mesh>& mesh,
               const std::vector<BoundaryCondition>& boundary_conditions);

    void Solve();
    void Integrate() const;

    // Setters
    void SetBoundaryConditions(
        const std::vector<BoundaryCondition>& boundary_conditions) {
        engine_->boundary_conditions = boundary_conditions;
    }
    void SetYoungsModulus(float E) const { engine_->modulus_of_elasticity = E; }
    void SetPoissonsRatio(float v) const { engine_->poissons_ratio = v; }
    void SetMass(float mass) const { integrator_->SetMassMatrix(mass); }

  private:
    std::unique_ptr<LinearTetrahedral> engine_;
    std::unique_ptr<ExplicitCentralDifferenceMethod> integrator_;
};
