#pragma once

#include "BoundaryCondition.h"
#include "ExplicitCentralDifference.h"
#include "LinearTetrahedral.h"
#include "Mesh.h"

#include <memory>
#include <vector>

class Simulation {
  public:
    float dt = 0.01f;
    float current_time = 0.f;

    Simulation() = default;
    Simulation(float youngs_modulus, float poissons_ratio, float point_mass,
               const std::shared_ptr<Mesh>& mesh,
               const std::vector<BoundaryCondition>& boundary_conditions);

    void Solve();
    Eigen::VectorXf Integrate() const;

    // Setters
    void SetBoundaryConditions(
        const std::vector<BoundaryCondition>& boundary_conditions) {
        engine_->boundary_conditions = boundary_conditions;
    }
    void SetYoungsModulus(float E) { youngs_modulus_ = E; }
    void SetPoissonsRatio(float v) { poissons_ratio_ = v; }
    void SetMass(float mass) const { integrator_->SetMassMatrix(mass); }
    void SetNodalMass(float mass) const { integrator_->SetMassMatrix(mass); }
    void SetTimestepSize(float dt) { this->dt = dt; }

    void SetRayleighMu(float mu) {
        rayleigh_mu_ = mu;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetRayleighLambda(float lambda) {
        rayleigh_lambda_ = lambda;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetDampingMatrix(float mu, float lambda) {
        rayleigh_mu_ = mu;
        rayleigh_lambda_ = lambda;
        integrator_->SetDamping(rayleigh_mu_, rayleigh_lambda_);
    }

    void SetTetgenFlags(const std::string& value) const {
        mesh_->SetTetgenFlags(value);
    }

    void SetSliceValue(float value) const { mesh_->SetCutPlane(value); }
    void SetSliceAxis(const std::string& value) const {
        mesh_->SetCutPlaneAxis(StringToCutPlaneAxis(value));
    }

    void SetMesh(const std::shared_ptr<Mesh>& mesh) { mesh_ = mesh; }

    // Getters
    [[nodiscard]] float SliceValue() const { return mesh_->cut_plane; }
    [[nodiscard]] float NodalMass() const { return integrator_->NodalMass(); }
    [[nodiscard]] float PoissonsRatio() const { return poissons_ratio_; }
    [[nodiscard]] float YoungsModulus() const { return youngs_modulus_; }
    [[nodiscard]] float TimestepSize() const { return dt; }

    [[nodiscard]] float RayleighLambda() const { return rayleigh_lambda_; }
    [[nodiscard]] float RayleighMu() const { return rayleigh_mu_; }

  private:
    float rayleigh_mu_ = 0.5f;
    float rayleigh_lambda_ = 0.5f;

    float youngs_modulus_ = 10000.f;
    float poissons_ratio_ = 0.3f;

    std::unique_ptr<LinearTetrahedral> engine_;
    std::unique_ptr<ExplicitCentralDifferenceMethod> integrator_;
    std::shared_ptr<Mesh> mesh_;
};
