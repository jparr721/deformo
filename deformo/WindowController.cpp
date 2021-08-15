#include "WindowController.h"
#include "Mesh.h"
#include "QTUtils.h"
#include "Utils.h"

WindowController::WindowController(Ui::deformoClass& ui) : ui_(ui) {

    const std::string default_mesh_path = "cube.ply";
    mesh = std::make_shared<Mesh>(default_mesh_path, tetgen_flags_);

    simulation_controller_ = std::make_unique<SimulationController>(mesh);
    designer_controller_ = std::make_unique<DesignerController>(mesh);

    ConnectUiElementsToSimulation();
    ConnectUiElementsToDesigner();
    ui_.designer_dataset_generator_progressbar->setVisible(false);
}

void WindowController::SetRenderer(const std::shared_ptr<Renderer>& renderer) {
    renderer_ = renderer;
}

void WindowController::StepForward() { simulation_controller_->StepForward(); }

void WindowController::Reset() { simulation_controller_->Reset(ui_); }

bool WindowController::IsSimulating() {
    return simulation_controller_->simulating;
}

void WindowController::SetSliceAxis(const QString& value) {
    slice_axis_ = utils::qt::QStringToString(value);
    emit OnSliceAxisChange(value);
}

void WindowController::SetSliceValue(Real value) {
    mesh->SetSliceValue(value);
    emit OnSliceValueChange(value);
}

void WindowController::SetNodalMass(Real value) {
    simulation_controller_->SetNodalMass(value);
    emit OnNodalMassChange(value);
}

void WindowController::SetPoissonsRatio(double value) {
    simulation_controller_->SetPoissonsRatio(value);
    emit OnPoissonsRatioChange(value);
}

void WindowController::SetYoungsModulus(double value) {
    simulation_controller_->SetYoungsModulus(value);
    emit OnYoungsModulusChange(value);
}

void WindowController::SetTimestepSize(double value) {
    simulation_controller_->SetTimestepSize(value);
    emit OnTimestepSizeChange(value);
}

void WindowController::SetRayleighLambda(double value) {
    simulation_controller_->SetRayleighLambda(value);
    emit OnRayleighLambdaChange(value);
}

void WindowController::SetRayleighMu(double value) {
    simulation_controller_->SetRayleighMu(value);
    emit OnRayleighMuChange(value);
}

void WindowController::RunSimulationButtonPressed() {
    simulation_controller_->RunSimulationButtonPressed(ui_);
}

void WindowController::SetTetgenFlags(const QString& value) {
    tetgen_flags_ = utils::qt::QStringToString(value);
    emit OnTetgenFlagsChange(value);
}

void WindowController::SetRenderMode(bool checked) {
    if (checked) {
        render_mode_ = GL_FILL;
    } else {
        render_mode_ = GL_LINE;
    }
}

void WindowController::RenderSimulationButtonPressed() {
    Reset();
    renderer_->SetRenderMode(render_mode_);

    mesh->SetTetgenFlags(tetgen_flags_);
}

void WindowController::PlaybackSkipStartButtonPressed() {
    simulation_controller_->PlaybackSkipStartButtonPressed(ui_);
}

void WindowController::PlaybackSkipEndButtonPressed() {
    simulation_controller_->PlaybackSkipEndButtonPressed(ui_);
}

void WindowController::PlaybackPauseButtonPressed() {
    simulation_controller_->PlaybackPauseButtonPressed(ui_);
}

void WindowController::PlaybackPlayButtonPressed() {
    simulation_controller_->PlaybackPlayButtonPressed();
}

void WindowController::PlaybackSliderChanged(int value) {
    simulation_controller_->PlaybackSliderChanged(value);
}

void WindowController::TabBarClicked(int index) {
    simulation_controller_->TabBarClicked(index, ui_);
}

void WindowController::SetForceSquareDimensions(bool checked) {
    if (checked) {
        ui_.designer_width_spinbox->setDisabled(true);
        ui_.designer_depth_spinbox->setDisabled(true);
        ui_.designer_height_spinbox_label->setText("Dimensions (X, Y, Z):");
    } else {
        ui_.designer_width_spinbox->setDisabled(false);
        ui_.designer_depth_spinbox->setDisabled(false);
        ui_.designer_height_spinbox_label->setText("Height:");
    }
}

void WindowController::SetImplicitSurfaceHeight(int value) {
    designer_controller_->SetImplicitSurfaceHeight(value);
    emit OnSetImplicitSurfaceHeight(value);
}

void WindowController::SetImplicitSurfaceWidth(int value) {
    designer_controller_->SetImplicitSurfaceWidth(value);
    emit OnSetImplicitSurfaceWidth(value);
}

void WindowController::SetImplicitSurfaceDepth(int value) {
    designer_controller_->SetImplicitSurfaceDepth(value);
    emit OnSetImplicitSurfaceDepth(value);
}

void WindowController::SetSquareShapedMaterial(bool checked) {
    designer_controller_->SetSquareShapedMaterial(checked);
}

void WindowController::SetUniformMaterial(bool checked) {
    designer_controller_->SetUniformMaterial(checked, ui_);
}

void WindowController::SetIsotropicMaterial(bool checked) {
    designer_controller_->SetIsotropicMaterial(checked);
}

void WindowController::SetNumberOfInclusions(int value) {
    designer_controller_->SetNumberOfInclusions(value);
    emit OnSetNumberOfInclusions(value);
}

void WindowController::SetMaterialOneName(const QString& value) {
    const std::string name = utils::qt::QStringToString(value);
    designer_controller_->SetMaterialOneName(name);
    emit OnSetMaterialOneName(value);
}

void WindowController::SetMaterialOnePoissonsRatio(double value) {
    Real v = value;
    designer_controller_->SetMaterialOnePoissonsRatio(value);
    emit OnSetMaterialOnePoissonsRatio(value);
}

void WindowController::SetMaterialOneYoungsModulus(double value) {
    Real E = value;
    designer_controller_->SetMaterialOneYoungsModulus(value);
    emit OnSetMaterialOneYoungsModulus(value);
}

void WindowController::SetMaterialTwoName(const QString& value) {
    const std::string name = utils::qt::QStringToString(value);
    designer_controller_->SetMaterialTwoName(name);
    emit OnSetMaterialTwoName(value);
}

void WindowController::SetMaterialTwoPoissonsRatio(double value) {
    Real v = value;
    designer_controller_->SetMaterialTwoPoissonsRatio(value);
    emit OnSetMaterialTwoPoissonsRatio(value);
}

void WindowController::SetMaterialTwoYoungsModulus(double value) {
    Real E = value;
    designer_controller_->SetMaterialTwoYoungsModulus(value);
    emit OnSetMaterialTwoYoungsModulus(value);
}

void WindowController::ComputeDesignedShapeButtonPressed() {
    designer_controller_->ComputeDesignedShapeButtonPressed();
}

void WindowController::SetInclusionHeight(int value) {
    designer_controller_->SetInclusionHeight(value);
    emit OnSetInclusionHeight(value);
}

void WindowController::SetInclusionWidth(int value) {
    designer_controller_->SetInclusionWidth(value);
    emit OnSetInclusionWidth(value);
}

void WindowController::SetInclusionDepth(int value) {
    designer_controller_->SetInclusionDepth(value);
    emit OnSetInclusionDepth(value);
}

void WindowController::SetSquareShapedInclusion(bool checked) {
    designer_controller_->SetSquareShapedInclusion(checked);
}

void WindowController::ConnectUiElementsToSimulation() {
    // Slice Axis Combo Box
    ui_.slice_axis_combo_box->addItem("X-Axis");
    ui_.slice_axis_combo_box->addItem("Y-Axis");
    ui_.slice_axis_combo_box->addItem("Z-Axis");
    connect(ui_.slice_axis_combo_box, &QComboBox::currentTextChanged, this,
            &WindowController::SetSliceAxis);
    connect(this, &WindowController::OnSliceAxisChange,
            ui_.slice_axis_combo_box, &QComboBox::setCurrentText);

    // Slice Axis Slider
    connect(ui_.slice_value_slider, &QSlider::valueChanged, this,
            &WindowController::SetSliceValue);
    connect(this, &WindowController::OnSliceValueChange, ui_.slice_value_slider,
            &QSlider::setValue);

    // Timestep Spin Box
    connect(ui_.timestep_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetTimestepSize);
    connect(this, &WindowController::OnTimestepSizeChange,
            ui_.timestep_double_spin_box, &QDoubleSpinBox::setValue);

    // Nodal Mass Spin Box
    connect(ui_.nodal_mass_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetNodalMass);
    connect(this, &WindowController::OnNodalMassChange,
            ui_.nodal_mass_double_spin_box, &QDoubleSpinBox::setValue);

    // Poissons Ratio Spin Box
    connect(ui_.poissons_ratio_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetPoissonsRatio);
    connect(this, &WindowController::OnPoissonsRatioChange,
            ui_.poissons_ratio_double_spin_box, &QDoubleSpinBox::setValue);

    // Young's Modulus Spin Box
    connect(ui_.youngs_modulus_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetYoungsModulus);
    connect(this, &WindowController::OnYoungsModulusChange,
            ui_.youngs_modulus_double_spin_box, &QDoubleSpinBox::setValue);

    // Damping Parameters
    // Rayleigh Damping Lambda
    connect(ui_.rayleigh_lambda_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetRayleighLambda);
    connect(this, &WindowController::OnRayleighLambdaChange,
            ui_.rayleigh_lambda_double_spin_box, &QDoubleSpinBox::setValue);

    // Rayleigh Damping Mu
    connect(ui_.rayleigh_mu_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetRayleighMu);
    connect(this, &WindowController::OnRayleighMuChange,
            ui_.rayleigh_mu_double_spin_box, &QDoubleSpinBox::setValue);

    // Run Simulation Button (Right Pane, Sim Settings)
    connect(ui_.sim_settings_run_button, &QPushButton::released, this,
            &WindowController::RunSimulationButtonPressed);

    // Render Parameters
    connect(ui_.tetgen_flags_line_edit, &QLineEdit::textChanged, this,
            &WindowController::SetTetgenFlags);
    connect(this, &WindowController::OnTetgenFlagsChange,
            ui_.tetgen_flags_line_edit, &QLineEdit::setText);

    connect(ui_.render_filled_radio_button, &QRadioButton::clicked, this,
            &WindowController::SetRenderMode);

    // Re Render Simulation (Right Pane, Render Settings)
    connect(ui_.render_properties_render_button, &QPushButton::released, this,
            &WindowController::RenderSimulationButtonPressed);

    // Playback Buttons
    connect(ui_.play_push_button, &QPushButton::released, this,
            &WindowController::PlaybackPlayButtonPressed);
    connect(ui_.pause_push_button, &QPushButton::released, this,
            &WindowController::PlaybackPauseButtonPressed);
    connect(ui_.skip_beginning_push_button, &QPushButton::released, this,
            &WindowController::PlaybackSkipStartButtonPressed);
    connect(ui_.skip_end_push_button, &QPushButton::released, this,
            &WindowController::PlaybackSkipEndButtonPressed);

    // Playback Slider
    connect(ui_.playback_controller, &QSlider::valueChanged, this,
            &WindowController::PlaybackSliderChanged);
    connect(this, &WindowController::OnPlaybackSliderChange,
            ui_.playback_controller, &QSlider::setValue);

    // Tab Bar Affects Sim State
    connect(ui_.tabWidget, &QTabWidget::tabBarClicked, this,
            &WindowController::TabBarClicked);

    // Initializers for UI Components
    SetSliceValue(simulation_controller_->simulation->SliceValue());
    SetNodalMass(simulation_controller_->simulation->NodalMass());
    SetPoissonsRatio(simulation_controller_->simulation->PoissonsRatio());
    SetYoungsModulus(simulation_controller_->simulation->YoungsModulus());
    SetTimestepSize(simulation_controller_->simulation->TimestepSize());
    SetTetgenFlags(QString::fromUtf8(tetgen_flags_.c_str()));

    SetRayleighLambda(simulation_controller_->simulation->RayleighLambda());
    SetRayleighMu(simulation_controller_->simulation->RayleighMu());
}

void WindowController::ConnectUiElementsToDesigner() {
    // First 3 checkboxes
    connect(ui_.designer_force_square_dimensions_checkbox, &QCheckBox::clicked,
            this, &WindowController::SetForceSquareDimensions);
    connect(ui_.designer_uniform_material_checkbox, &QCheckBox::clicked, this,
            &WindowController::SetUniformMaterial);
    connect(ui_.designer_isotropic_material_checkbox, &QCheckBox::clicked, this,
            &WindowController::SetIsotropicMaterial);

    // Implicit Surface Values
    connect(ui_.designer_height_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetImplicitSurfaceHeight);
    connect(this, &WindowController::OnSetImplicitSurfaceHeight,
            ui_.designer_height_spinbox, &QSpinBox::setValue);

    connect(ui_.designer_width_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetImplicitSurfaceWidth);
    connect(this, &WindowController::OnSetImplicitSurfaceWidth,
            ui_.designer_width_spinbox, &QSpinBox::setValue);

    connect(ui_.designer_depth_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetImplicitSurfaceDepth);
    connect(this, &WindowController::OnSetImplicitSurfaceDepth,
            ui_.designer_depth_spinbox, &QSpinBox::setValue);

    // Inclusions Menu Section
    connect(ui_.designer_material_2_number_of_inclusions,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetNumberOfInclusions);
    connect(this, &WindowController::OnSetNumberOfInclusions,
            ui_.designer_material_2_number_of_inclusions, &QSpinBox::setValue);

    connect(ui_.designer_inclusion_force_square_dimensions_checkbox,
            &QCheckBox::clicked, this,
            &WindowController::SetSquareShapedInclusion);

    connect(ui_.designer_inclusion_height_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetInclusionHeight);
    connect(this, &WindowController::OnSetInclusionHeight,
            ui_.designer_inclusion_height_spinbox, &QSpinBox::setValue);

    connect(ui_.designer_inclusion_width_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetInclusionWidth);
    connect(this, &WindowController::OnSetInclusionWidth,
            ui_.designer_inclusion_width_spinbox, &QSpinBox::setValue);

    connect(ui_.designer_inclusion_depth_spinbox,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &WindowController::SetInclusionDepth);
    connect(this, &WindowController::OnSetInclusionDepth,
            ui_.designer_inclusion_depth_spinbox, &QSpinBox::setValue);

    // Materials Menu Section
    connect(ui_.designer_material_1_name_line_edit, &QLineEdit::textChanged,
            this, &WindowController::SetMaterialOneName);
    connect(this, &WindowController::OnSetMaterialOneName,
            ui_.designer_material_1_name_line_edit, &QLineEdit::setText);

    connect(ui_.designer_material_1_poissons_ratio_double_spinbox,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetMaterialOnePoissonsRatio);
    connect(this, &WindowController::OnSetMaterialOnePoissonsRatio,
            ui_.designer_material_1_poissons_ratio_double_spinbox,
            &QDoubleSpinBox::setValue);

    connect(ui_.designer_material_1_youngs_modulus_double_spinbox,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetMaterialOneYoungsModulus);
    connect(this, &WindowController::OnSetMaterialOneYoungsModulus,
            ui_.designer_material_1_youngs_modulus_double_spinbox,
            &QDoubleSpinBox::setValue);

    connect(ui_.designer_material_2_name_line_edit, &QLineEdit::textChanged,
            this, &WindowController::SetMaterialTwoName);
    connect(this, &WindowController::OnSetMaterialTwoName,
            ui_.designer_material_2_name_line_edit, &QLineEdit::setText);

    connect(ui_.designer_material_2_poissons_ratio_double_spinbox,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetMaterialTwoPoissonsRatio);
    connect(this, &WindowController::OnSetMaterialTwoPoissonsRatio,
            ui_.designer_material_2_poissons_ratio_double_spinbox,
            &QDoubleSpinBox::setValue);

    connect(ui_.designer_material_2_youngs_modulus_double_spinbox,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetMaterialTwoYoungsModulus);
    connect(this, &WindowController::OnSetMaterialTwoYoungsModulus,
            ui_.designer_material_2_youngs_modulus_double_spinbox,
            &QDoubleSpinBox::setValue);

    // Compute Button
    connect(ui_.designer_compute_shape_push_button, &QPushButton::pressed, this,
            &WindowController::ComputeDesignedShapeButtonPressed);
}

void WindowController::ResetPlaybackControls() {
    simulation_controller_->ResetPlaybackControls(ui_);
}
