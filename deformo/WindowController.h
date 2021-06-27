#pragma once

#include "BoundaryCondition.h"
#include "Simulation.h"
#include "ui_deformo.h"
#include <QObject>
#include <QString>
#include <memory>
#include <string>
#include <vector>

class WindowController : public QObject {
    Q_OBJECT
  public:
    unsigned int steps_taken = 0;
    unsigned int max_steps = 0;

    std::shared_ptr<Mesh> mesh;

    BoundaryConditions boundary_conditions;

    WindowController(Ui::deformoClass& ui, const std::string& mesh_path);

    void StepForward();
    void StepBackward();

    void Reset();

  public slots:
    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void SetSliceAxis(const QString& value);
    void SetSliceValue(float value);
    void SetNodalMass(float value);
    void SetPoissonsRatio(double value);
    void SetYoungsModulus(double value);
    void SetTimestepSize(double value);

    // Simulation Settings -- Damping Parameters
    void SetRayleighLambda(double value);
    void SetRayleighMu(double value);

    // Simulation Settings -- Run Simulation Button
    void RunSimulationButtonPressed();

    // Render Properties Window -- Tetgen Flags
    void SetTetgenFlags(const QString& value);

    // Render Properties Window -- Re Render Button
    void RenderSimulationButtonPressed();

  signals:
    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void OnSliceAxisChange(const QString& value);
    void OnSliceValueChange(float value);
    void OnNodalMassChange(float value);
    void OnPoissonsRatioChange(double value);
    void OnYoungsModulusChange(double value);
    void OnTimestepSizeChange(double value);

    // Simulation Settings -- Damping Parameters
    void OnRayleighLambdaChange(double value);
    void OnRayleighMuChange(double value);

    // Render Properties Window
    void OnTetgenFlagsChange(const QString& value);

  private:
    // Simulation Runtime Parameters
    float dt_ = 0.01f;

    // Physical Parameters
    float nodal_mass_ = 1.f;

    // Damping Parameters
    float rayleigh_mu_ = 0.5f;
    float rayleigh_lambda_ = 0.5f;

    // FEA Parameters
    float youngs_modulus_ = 10000.f;
    float poissons_ratio_ = 0.3f;

    // Render Parameters
    std::string tetgen_flags_;

    // Render-Mutation Parameters
    float slice_value_;
    std::string slice_axis_;

    std::unique_ptr<Simulation> simulation_;

    BoundaryConditions
    GenerateDefaultBoundaryConditions(const std::shared_ptr<Mesh>& mesh);

  private:
    void ConnectUiElementsToSimulation(Ui::deformoClass& ui);
};
