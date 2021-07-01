#pragma once

#include "BoundaryCondition.h"
#include "ExplicitCentralDifference.h"
#include "LinearTetrahedral.h"
#include "Mesh.h"

#include <memory>
#include <vector>

class Simulation {
  public:
    Real dt = 0.01f;
    Real current_time = 0.f;

    Simulation() = default;
    Simulation(Real youngs_modulus, Real poissons_ratio, Real point_mass,
               const std::shared_ptr<Mesh>& mesh,
               const std::vector<BoundaryCondition>& boundary_conditions);

    void Solve();
    VectorXr Integrate() const;

    // Setters
    void SetBoundaryConditions(
        const std::vector<BoundaryCondition>& boundary_conditions) {
        engine_->boundary_conditions = boundary_conditions;
    }
    void SetYoungsModulus(Real E) { youngs_modulus_ = E; }
    void SetPoissonsRatio(Real v) { poissons_ratio_ = v; }
    void SetMass(Real mass) const { integrator_->SetMassMatrix(mass); }
    void SetNodalMass(Real mass) const { integrator_->SetMassMatrix(mass); }
    void SetTimestepSize(Real dt) { this->dt = dt; }

    void SetRayleighMu(Real mu) {
        rayleigh_mu_ = mu;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetRayleighLambda(Real lambda) {
        rayleigh_lambda_ = lambda;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetDampingMatrix(Real mu, Real lambda) {
        rayleigh_mu_ = mu;
        rayleigh_lambda_ = lambda;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetTetgenFlags(const std::string& value) const {
        mesh_->SetTetgenFlags(value);
    }

    void SetSliceValue(Real value) const { mesh_->SetSliceValue(value); }
    void SetSliceAxis(const std::string& value) const {
        mesh_->SetSliceAxis(StringToSliceAxis(value));
    }

    void SetMesh(const std::shared_ptr<Mesh>& mesh) { mesh_ = mesh; }

    // Getters
    [[nodiscard]] Real SliceValue() const { return mesh_->slice_value; }
    [[nodiscard]] Real NodalMass() const { return integrator_->NodalMass(); }
    [[nodiscard]] Real PoissonsRatio() const { return poissons_ratio_; }
    [[nodiscard]] Real YoungsModulus() const { return youngs_modulus_; }
    [[nodiscard]] Real TimestepSize() const { return dt; }

    [[nodiscard]] Real RayleighLambda() const { return rayleigh_lambda_; }
    [[nodiscard]] Real RayleighMu() const { return rayleigh_mu_; }

  private:
    Real rayleigh_mu_ = 0.5f;
    Real rayleigh_lambda_ = 0.5f;

    Real youngs_modulus_ = 10000.f;
    Real poissons_ratio_ = 0.3f;

    std::unique_ptr<LinearTetrahedral> engine_;
    std::unique_ptr<ExplicitCentralDifferenceMethod> integrator_;
    std::shared_ptr<Mesh> mesh_;
};
