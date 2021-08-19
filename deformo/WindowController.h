#pragma once

#include "BoundaryCondition.h"
#include "DesignerController.h"
#include "Simulation.h"
#include "SimulationController.h"
#include "ui_deformo.h"
#include <QObject>
#include <QString>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class WindowController : public QObject {
    Q_OBJECT
  public:
    unsigned int steps_taken = 0;
    unsigned int max_steps = 0;

    std::shared_ptr<Mesh> mesh;

    BoundaryConditions boundary_conditions;

    WindowController(Ui::deformoClass& ui);

    void SetRenderer(const std::shared_ptr<Renderer>& renderer);

    void StepForward();

    void Reset();

    bool IsSimulating();

  public slots:
    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void SetSliceAxis(const QString& value);
    void SetSliceValue(Real value);
    void SetNodalMass(Real value);
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

    // ====== Designer Controls Window
    void TabBarClicked(int index);

    // Designer -- Dimensions
    void SetForceSquareDimensions(bool checked);

    // Designer -- Size of shape
    void SetImplicitSurfaceHeight(int value); // Rows
    void SetImplicitSurfaceWidth(int value);  // Cols
    void SetImplicitSurfaceDepth(int value);  // Layers

    // Designer -- Square Material
    void SetSquareShapedMaterial(bool checked);

    // Designer -- 1-Material
    void SetUniformMaterial(bool checked);

    // Designer -- Isotropic Material
    void SetIsotropicMaterial(bool checked);

    // Designer -- Size of inclusion
    void SetNumberOfInclusions(int value);

    // Designer -- Material 1 Specifications
    void SetMaterialOneName(const QString& value);
    void SetMaterialOnePoissonsRatio(double value);
    void SetMaterialOneYoungsModulus(double value);

    // Designer -- Material 2 Specifications
    void SetMaterialTwoName(const QString& value);
    void SetMaterialTwoPoissonsRatio(double value);
    void SetMaterialTwoYoungsModulus(double value);

    // Designer -- Compute Button
    void ComputeDesignedShapeButtonPressed();

    // Designer -- Inclusion Shape Options
    void SetInclusionHeight(int value);
    void SetInclusionWidth(int value);
    void SetInclusionDepth(int value);

    // Designer -- Inclusion Is Square.
    void SetSquareShapedInclusion(bool checked);

    // Designer -- Dataset Generator -- CSV Path
    void SetCSVPathButtonClicked();
    void SetOutputCSVFileName(const QString& value);
    void SetGenerator(const QString& value);
    void SetDatasetGeneratorNumberOfEntries(int value);
    void DatasetGeneratorComputeButtonPressed();


  signals:
    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void OnSliceAxisChange(const QString& value);
    void OnSliceValueChange(Real value);
    void OnNodalMassChange(Real value);
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

    // ====== Designer Controls Window
    // Implicit surface
    void OnSetImplicitSurfaceHeight(int value);
    void OnSetImplicitSurfaceWidth(int value);
    void OnSetImplicitSurfaceDepth(int value);

    // Designer -- Size of inclusion
    void OnSetNumberOfInclusions(int value);

    // Designer -- Material 1 Specifications
    void OnSetMaterialOneName(const QString& value);
    void OnSetMaterialOnePoissonsRatio(double value);
    void OnSetMaterialOneYoungsModulus(double value);

    // Designer -- Material 2 Specifications
    void OnSetMaterialTwoName(const QString& value);
    void OnSetMaterialTwoPoissonsRatio(double value);
    void OnSetMaterialTwoYoungsModulus(double value);

    // Designer -- Inclusion Shape Options
    void OnSetInclusionHeight(int value);
    void OnSetInclusionWidth(int value);
    void OnSetInclusionDepth(int value);

    // Designer -- Dataset Generatr -- CSV Path
    void OnSetOutputCSVPath(const QString& value);
    void OnSetOutputCSVFileName(const QString& value);
    void OnSetGenerator(const QString& value);

  private:
    Ui::deformoClass ui_;

    // =========== Simulator ==============
    std::unique_ptr<SimulationController> simulation_controller_;
    std::unique_ptr<DesignerController> designer_controller_;

    GLenum render_mode_ = GL_LINES;

    // Render Parameters
    std::string tetgen_flags_ = "zpq";

    // Render-Mutation Parameters
    Real slice_value_;
    std::string slice_axis_;

    std::shared_ptr<Renderer> renderer_;

    // =========== Simulator ==============

    void ConnectUiElementsToSimulation();
    void ConnectUiElementsToDesigner();
    void ResetPlaybackControls();
};
