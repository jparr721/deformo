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

    void SetRenderer(const std::shared_ptr<Renderer> renderer);

    void StepForward();

    void Reset();

    bool IsSimulating();

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

    // Render Properties Window -- Render Mode
    void SetRenderMode(bool checked);

    // Render Properties Window -- Re Render Button
    void RenderSimulationButtonPressed();

    // Playback Controls
    void PlaybackSkipStartButtonPressed();
    void PlaybackSkipEndButtonPressed();
    void PlaybackPauseButtonPressed();
    void PlaybackPlayButtonPressed();
    void PlaybackSliderChanged(int value);

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

    // Playback Controls
    void OnPlaybackSliderChange(int value);

  private:
    bool simulating_ = false;
    // Simulation Runtime Parameters
    float dt_ = 0.01f;

    // Physical Parameters
    float nodal_mass_ = 5.f;

    // Damping Parameters
    float rayleigh_mu_ = 0.5f;
    float rayleigh_lambda_ = 0.5f;

    // FEA Parameters
    float youngs_modulus_ = 10000.f;
    float poissons_ratio_ = 0.3f;

    GLenum render_mode_ = GL_LINES;

    // Render Parameters
    std::string tetgen_flags_ = "zpq";
    // Render-Mutation Parameters
    float slice_value_;
    std::string slice_axis_;

    std::unique_ptr<Simulation> simulation_;

    std::shared_ptr<Renderer> renderer_;

    std::vector<Eigen::VectorXf> recorded_displacements_;

    BoundaryConditions
    GenerateDefaultBoundaryConditions(const std::shared_ptr<Mesh>& mesh);

    Ui::deformoClass ui_;

  private:
    void ConnectUiElementsToSimulation();
    void DisableStaticUiElements();
    void EnableStaticUiElements();
    void ResetPlaybackControls();
    void RecomputeSliceValueRange();
};
