#pragma once

#include "BoundaryCondition.h"
#include "Simulation.h"
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

    // Designer -- 1-Material
    void SetUniformMaterial(bool checked);

    // Designer -- Isotropic Material
    void SetIsotropicMaterial(bool checked);

    // Designer -- Size of inclusion
    void SetNumberOfInclusions(int value);

    // Designer -- Material 1 Specifications
    void SetMaterialOneName(const QString& value);
    void SetMaterialOnePoissionsRatio(double value);
    void SetMaterialOneYoungsModulus(double value);

    // Designer -- Material 2 Specifications
    void SetMaterialTwoName(const QString& value);
    void SetMaterialTwoPoissionsRatio(double value);
    void SetMaterialTwoYoungsModulus(double value);

    // Designer -- Compute Button
    void ComputeDesignedShapeButtonPressed();

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
    void OnSetMaterialOnePoissionsRatio(double value);
    void OnSetMaterialOneYoungsModulus(double value);

    // Designer -- Material 2 Specifications
    void OnSetMaterialTwoName(const QString& value);
    void OnSetMaterialTwoPoissionsRatio(double value);
    void OnSetMaterialTwoYoungsModulus(double value);

  private:
    enum TabWindow {
        kDesigner = 0x00,
        kSimulation = 0x01,
    };

    Ui::deformoClass ui_;

    // =========== Simulator ==============
    bool simulating_ = false;
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

    GLenum render_mode_ = GL_LINES;

    // Render Parameters
    std::string tetgen_flags_ = "zpq";
    // Render-Mutation Parameters
    Real slice_value_;
    std::string slice_axis_;

    std::unique_ptr<Simulation> simulation_;

    std::shared_ptr<Renderer> renderer_;

    std::vector<VectorXr> recorded_displacements_;

    BoundaryConditions
    GenerateDefaultBoundaryConditions(const std::shared_ptr<Mesh>& mesh);

    // =========== Simulator ==============

    // =========== Designer ==============

    // Implicit Surface Options
    int designer_implicit_surface_height = 0;
    int designer_implicit_surface_width = 0;
    int designer_implicit_surface_depth = 0;

    // Is material made of only material 1
    bool designer_is_uniform = false;

    // Isotropic Material Generator
    bool designer_is_isotropic = false;

    // Inclusion ratio for material 2
    int designer_material_2_number_of_inclusions = 0;

    // Material 1
    Material material_1;

    // Material 2
    Material material_2;

    // =========== Designer ==============

    void ConnectUiElementsToSimulation();
    void ConnectUiElementsToDesigner();
    void DisableStaticUiElements();
    void EnableStaticUiElements();
    void ResetPlaybackControls();
    void RecomputeSliceValueRange();
};
