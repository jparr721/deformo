#include "Window.h"

Window::Window(QWidget* parent) : QMainWindow(parent) {
    ui.setupUi(this);
    window_controller_ = new WindowController(ui, "cube.ply");
    connect(ui.sim_settings_run_button, &QPushButton::released,
            window_controller_, &WindowController::RunSimulationButtonPressed);

    // Slice Axis Combo Box
    ui.slice_axis_combo_box->addItem("X-Axis");
    ui.slice_axis_combo_box->addItem("Y-Axis");
    ui.slice_axis_combo_box->addItem("Z-Axis");
    connect(ui.slice_axis_combo_box, &QComboBox::currentTextChanged,
            ui.openGLWidget, &GLWidget::SetSliceAxis);
    connect(ui.openGLWidget, &GLWidget::OnSliceAxisChange,
            ui.slice_axis_combo_box, &QComboBox::setCurrentText);

    // Slice Axis Slider
    connect(ui.slice_value_slider, &QSlider::valueChanged, ui.openGLWidget,
            &GLWidget::SetSliceValue);
    connect(ui.openGLWidget, &GLWidget::OnSliceValueChange,
            ui.slice_value_slider, &QSlider::setValue);

    // Timestep Spin Box
    connect(ui.timestep_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetTimestepSize);
    connect(ui.openGLWidget, &GLWidget::OnTimestepSizeChange,
            ui.timestep_double_spin_box, &QDoubleSpinBox::setValue);

    // Nodal Mass Spin Box
    connect(ui.nodal_mass_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetNodalMass);
    connect(ui.openGLWidget, &GLWidget::OnNodalMassChange,
            ui.nodal_mass_double_spin_box, &QDoubleSpinBox::setValue);

    // Poissons Ratio Spin Box
    connect(ui.poissons_ratio_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetPoissonsRatio);
    connect(ui.openGLWidget, &GLWidget::OnPoissonsRatioChange,
            ui.poissons_ratio_double_spin_box, &QDoubleSpinBox::setValue);

    // Young's Modulus Spin Box
    connect(ui.youngs_modulus_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetYoungsModulus);
    connect(ui.openGLWidget, &GLWidget::OnYoungsModulusChange,
            ui.youngs_modulus_double_spin_box, &QDoubleSpinBox::setValue);

    // Damping Parameters
    // Rayleigh Damping Lambda
    connect(ui.rayleigh_lambda_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetRayleighLambda);
    connect(ui.openGLWidget, &GLWidget::OnRayleighLambdaChange,
            ui.rayleigh_lambda_double_spin_box, &QDoubleSpinBox::setValue);

    // Rayleigh Damping Mu
    connect(ui.rayleigh_mu_double_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetRayleighMu);
    connect(ui.openGLWidget, &GLWidget::OnRayleighMuChange,
            ui.rayleigh_mu_double_spin_box, &QDoubleSpinBox::setValue);

    // Run Simulation Button (Right Pane, Sim Settings)
    // connect(ui.sim_settings_run_button, &QPushButton::released,
    // ui.openGLWidget,
    //        &GLWidget::RunSimulationButtonPressed);

    // Render Parameters
    connect(ui.tetgen_flags_line_edit, &QLineEdit::textChanged, ui.openGLWidget,
            &GLWidget::SetTetgenFlags);
    connect(ui.openGLWidget, &GLWidget::OnTetgenFlagsChange,
            ui.tetgen_flags_line_edit, &QLineEdit::setText);

    // Re Render Simulation (Right Pane, Render Settings)
    connect(ui.render_properties_render_button, &QPushButton::released,
            ui.openGLWidget, &GLWidget::RenderSimulationButtonPressed);
}

Window::~Window() { delete window_controller_; }
