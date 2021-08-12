#pragma once

#include "BoundaryCondition.h"
#include "Numerics.h"
#include "Simulation.h"
#include "ui_deformo.h"
#include <memory>
#include <string>

class SimulationController {
  public:
    bool simulating = false;

    std::unique_ptr<Simulation> simulation;

    SimulationController(const std::shared_ptr<Mesh>& mesh);

    void StepForward();

    void Reset(const Ui::deformoClass& ui);
    void ResetPlaybackControls(const Ui::deformoClass& ui);

    // Simulation Settings -- FEA Parameters
    void SetNodalMass(Real value);
    void SetPoissonsRatio(double value);
    void SetYoungsModulus(double value);
    void SetTimestepSize(double value);

    // Simulation Settings -- Damping Parameters
    void SetRayleighLambda(double value);
    void SetRayleighMu(double value);

    // Simulation Settings -- Run Simulation Button
    void RunSimulationButtonPressed(const Ui::deformoClass& ui);

    // Playback Controls
    void PlaybackSkipStartButtonPressed(const Ui::deformoClass& ui);
    void PlaybackSkipEndButtonPressed(const Ui::deformoClass& ui);
    void PlaybackPauseButtonPressed(const Ui::deformoClass& ui);
    void PlaybackPlayButtonPressed();
    void PlaybackSliderChanged(int value);

    void TabBarClicked(int index, const Ui::deformoClass& ui);

  private:
    enum TabWindow {
        kDesigner = 0x00,
        kSimulation = 0x01,
    };
    // Sim
    // Simulation Runtime Parameters
    Real dt_ = 0.01f;

    // Physical Parameters
    Real nodal_mass_ = 5.f;

    // Damping Parameters
    Real rayleigh_mu_ = 0.5f;
    Real rayleigh_lambda_ = 0.5f;

    // FEA Parameters
    Real youngs_modulus_ = 10000.f;
    Real poissons_ratio_ = 0.3f;

    std::vector<VectorXr> recorded_displacements_;

    std::shared_ptr<Mesh> mesh_;

    BoundaryConditions
    GenerateDefaultBoundaryConditions(const std::shared_ptr<Mesh>& mesh);
};
