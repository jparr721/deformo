#include "WindowController.h"
#include "Mesh.h"
#include "QTUtils.h"
#include "Utils.h"

WindowController::WindowController(Ui::deformoClass& ui) : ui_(ui) {

    const std::string default_mesh_path = "cube.ply";
    mesh = std::make_shared<Mesh>(default_mesh_path, tetgen_flags_);
    boundary_conditions = GenerateDefaultBoundaryConditions(mesh);
    simulation_ =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh, boundary_conditions);

    ConnectUiElementsToSimulation();
}

void WindowController::SetRenderer(const std::shared_ptr<Renderer>& renderer) {
    renderer_ = renderer;
}

bool WindowController::IsSimulating() { return simulating_; }

void WindowController::StepForward() {
    simulation_->Solve();
    recorded_displacements_.push_back(simulation_->Integrate());
    RecomputeSliceValueRange();
}

void WindowController::Reset() {
    mesh->Reset();
    simulation_ =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh, boundary_conditions);
    simulation_->dt = dt_;
    simulation_->current_time = 0.f;
    steps_taken = 0;
    max_steps = 0;
    ResetPlaybackControls();
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
    nodal_mass_ = value;
    emit OnNodalMassChange(value);
}

void WindowController::SetPoissonsRatio(double value) {
    poissons_ratio_ = value;
    emit OnPoissonsRatioChange(value);
}

void WindowController::SetYoungsModulus(double value) {
    youngs_modulus_ = value;
    emit OnYoungsModulusChange(value);
}

void WindowController::SetTimestepSize(double value) {
    dt_ = value;
    emit OnTimestepSizeChange(value);
}

void WindowController::SetRayleighLambda(double value) {
    rayleigh_lambda_ = value;
    emit OnRayleighLambdaChange(value);
}

void WindowController::SetRayleighMu(double value) {
    rayleigh_mu_ = value;
    emit OnRayleighMuChange(value);
}

void WindowController::RunSimulationButtonPressed() { Reset(); }

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

    // TODO(@jparr721) Re-generate new mesh and sim (more work than I originally
    // thought).
    // mesh->SetTetgenFlags(tetgen_flags_);
}

void WindowController::PlaybackSkipStartButtonPressed() {
    ui_.playback_controller->setValue(0);
}

void WindowController::PlaybackSkipEndButtonPressed() {
    if (!recorded_displacements_.empty()) {
        ui_.playback_controller->setValue(recorded_displacements_.size() - 1);
    }
}

void WindowController::PlaybackPauseButtonPressed() {
    EnableStaticUiElements();
    if (!recorded_displacements_.empty()) {
        ui_.playback_controller->setMaximum(recorded_displacements_.size() - 1);
    }
    simulating_ = false;
}

void WindowController::PlaybackPlayButtonPressed() {
    DisableStaticUiElements();
    simulating_ = true;
}

void WindowController::PlaybackSliderChanged(int value) {
    mesh->Update(recorded_displacements_.at(value));
}

BoundaryConditions WindowController::GenerateDefaultBoundaryConditions(
    const std::shared_ptr<Mesh>& mesh) {
    const Eigen::Vector3f force(0.f, -100.f, 0.f);

    std::vector<unsigned int> indices;
    utils::FindMaxVertices(indices, mesh->positions);
    return AssignBoundaryConditionToFixedNodes(indices, force);
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

    // Initializers for UI Components
    SetSliceValue(simulation_->SliceValue());
    SetNodalMass(simulation_->NodalMass());
    SetPoissonsRatio(simulation_->PoissonsRatio());
    SetYoungsModulus(simulation_->YoungsModulus());
    SetTimestepSize(simulation_->TimestepSize());
    SetTetgenFlags(QString::fromUtf8(tetgen_flags_.c_str()));

    SetRayleighLambda(simulation_->RayleighLambda());
    SetRayleighMu(simulation_->RayleighMu());
}

void WindowController::DisableStaticUiElements() {}

void WindowController::EnableStaticUiElements() {}

void WindowController::ResetPlaybackControls() {
    recorded_displacements_.clear();
    recorded_displacements_.resize(0);
    PlaybackPauseButtonPressed();
}

void WindowController::RecomputeSliceValueRange() {
    unsigned int i = 0;

    if (mesh->slice_axis == SliceAxis::x_axis) {
        i = 0;
    } else if (mesh->slice_axis == SliceAxis::y_axis) {
        i = 1;
    } else if (mesh->slice_axis == SliceAxis::z_axis) {
        i = 2;
    }

    Real minimum = 0.f;
    Real maximum = 0.f;

    for (int i = 0; i < mesh->positions.size(); ++i) {
        minimum = std::fmin(minimum, mesh->positions(i));
        maximum = std::fmax(maximum, mesh->positions(i));
    }

    ui_.slice_value_slider->setMinimum(minimum);
    ui_.slice_value_slider->setMaximum(maximum);
}
