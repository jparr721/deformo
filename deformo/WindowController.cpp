#include "WindowController.h"
#include "Mesh.h"
#include "QTUtils.h"
#include "Utils.h"

WindowController::WindowController(Ui::deformoClass& ui,
                                   const std::string& mesh_path) {
    const std::string suffix = ".ply";
    // TODO (@jparr721) - Make this a ui error
    assert(0 == mesh_path.compare(mesh_path.size() - suffix.size(),
                                  suffix.size(), suffix) &&
           "YOU CAN ONLY LOAD .PLY FILES");

    mesh = std::make_shared<Mesh>(mesh_path, 0.f);
    boundary_conditions = GenerateDefaultBoundaryConditions(mesh);
    simulation_ =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh, boundary_conditions);

    ConnectUiElementsToSimulation(ui);
}

void WindowController::StepForward() {
    simulation_->Solve();
    simulation_->Integrate();
}

void WindowController::StepBackward() {}

void WindowController::Reset() {
    mesh->Reset();
    simulation_ =
        std::make_unique<Simulation>(youngs_modulus_, poissons_ratio_,
                                     nodal_mass_, mesh, boundary_conditions);
    steps_taken = 0;
    max_steps = 0;
}

void WindowController::SetSliceAxis(const QString& value) {
    slice_axis_ = utils::qt::QStringToString(value);
    emit OnSliceAxisChange(value);
}

void WindowController::SetSliceValue(float value) {
    slice_value_ = value;
    emit OnSliceValueChange(value);
}

void WindowController::SetNodalMass(float value) {
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

void WindowController::RunSimulationButtonPressed() {
    std::cout << "Run Simulation Pressed" << std::endl;
}

void WindowController::SetTetgenFlags(const QString& value) {
    tetgen_flags_ = utils::qt::QStringToString(value);
    emit OnTetgenFlagsChange(value);
}

void WindowController::RenderSimulationButtonPressed() {
    std::cout << "Render Button Pressed" << std::endl;
}

BoundaryConditions WindowController::GenerateDefaultBoundaryConditions(
    const std::shared_ptr<Mesh>& mesh) {
    const Eigen::Vector3f force(0.f, -100.f, 0.f);

    std::vector<unsigned int> indices;
    utils::FindMaxVertices(indices, mesh->positions);
    return AssignBoundaryConditionToFixedNodes(indices, force);
}

void WindowController::ConnectUiElementsToSimulation(Ui::deformoClass& ui) {
    // Slice Axis Combo Box
    ui.slice_axis_combo_box->addItem("X-Axis");
    ui.slice_axis_combo_box->addItem("Y-Axis");
    ui.slice_axis_combo_box->addItem("Z-Axis");
    connect(ui.slice_axis_combo_box, &QComboBox::currentTextChanged, this,
            &WindowController::SetSliceAxis);
    connect(this, &WindowController::OnSliceAxisChange, ui.slice_axis_combo_box,
            &QComboBox::setCurrentText);

    // Slice Axis Slider
    connect(ui.slice_value_slider, &QSlider::valueChanged, this,
            &WindowController::SetSliceValue);
    connect(this, &WindowController::OnSliceValueChange, ui.slice_value_slider,
            &QSlider::setValue);

    // Timestep Spin Box
    connect(ui.timestep_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetTimestepSize);
    connect(this, &WindowController::OnTimestepSizeChange,
            ui.timestep_double_spin_box, &QDoubleSpinBox::setValue);

    // Nodal Mass Spin Box
    connect(ui.nodal_mass_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetNodalMass);
    connect(this, &WindowController::OnNodalMassChange,
            ui.nodal_mass_double_spin_box, &QDoubleSpinBox::setValue);

    // Poissons Ratio Spin Box
    connect(ui.poissons_ratio_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetPoissonsRatio);
    connect(this, &WindowController::OnPoissonsRatioChange,
            ui.poissons_ratio_double_spin_box, &QDoubleSpinBox::setValue);

    // Young's Modulus Spin Box
    connect(ui.youngs_modulus_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetYoungsModulus);
    connect(this, &WindowController::OnYoungsModulusChange,
            ui.youngs_modulus_double_spin_box, &QDoubleSpinBox::setValue);

    // Damping Parameters
    // Rayleigh Damping Lambda
    connect(ui.rayleigh_lambda_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetRayleighLambda);
    connect(this, &WindowController::OnRayleighLambdaChange,
            ui.rayleigh_lambda_double_spin_box, &QDoubleSpinBox::setValue);

    // Rayleigh Damping Mu
    connect(ui.rayleigh_mu_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            &WindowController::SetRayleighMu);
    connect(this, &WindowController::OnRayleighMuChange,
            ui.rayleigh_mu_double_spin_box, &QDoubleSpinBox::setValue);

    // Run Simulation Button (Right Pane, Sim Settings)
    connect(ui.sim_settings_run_button, &QPushButton::released, this,
            &WindowController::RunSimulationButtonPressed);

    // Render Parameters
    connect(ui.tetgen_flags_line_edit, &QLineEdit::textChanged, this,
            &WindowController::SetTetgenFlags);
    connect(this, &WindowController::OnTetgenFlagsChange,
            ui.tetgen_flags_line_edit, &QLineEdit::setText);

    // Re Render Simulation (Right Pane, Render Settings)
    connect(ui.render_properties_render_button, &QPushButton::released, this,
            &WindowController::RenderSimulationButtonPressed);

    SetSliceValue(simulation_->SliceValue());
    SetNodalMass(simulation_->NodalMass());
    SetPoissonsRatio(simulation_->PoissonsRatio());
    SetYoungsModulus(simulation_->YoungsModulus());
    SetTimestepSize(simulation_->TimestepSize());
    SetTetgenFlags(QString::fromUtf8(simulation_->TetgenFlags().c_str()));

    SetRayleighLambda(simulation_->RayleighLambda());
    SetRayleighMu(simulation_->RayleighMu());
}
