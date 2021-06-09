#include "Window.h"

Window::Window(QWidget* parent) : QMainWindow(parent) {
    ui.setupUi(this);

    /* Sliders */
    connect(ui.cut_plane_slider, &QSlider::valueChanged, ui.openGLWidget,
            &GLWidget::SetCutPlane);
    connect(ui.openGLWidget, &GLWidget::OnCutPlaneChange, ui.cut_plane_slider,
            &QSlider::setValue);

    connect(ui.tetgen_flags_line_edit, &QLineEdit::textChanged, ui.openGLWidget,
            &GLWidget::OnTetgenFlagsChange);
    connect(ui.openGLWidget, &GLWidget::OnTetgenFlagsChange,
            ui.tetgen_flags_line_edit, &QLineEdit::setText);

    /* Check Boxes */
    connect(ui.interactive_mode_toggle, &QCheckBox::toggled, ui.openGLWidget,
            &GLWidget::SetInteractiveModeToggle);
    connect(ui.openGLWidget, &GLWidget::OnInteractiveModeToggled,
            ui.interactive_mode_toggle, &QCheckBox::setEnabled);

    /* Spin Boxes */
    connect(ui.poissons_ratio_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetPoissonsRatio);
    connect(ui.openGLWidget, &GLWidget::OnPoissonsRatioChange,
            ui.poissons_ratio_spin_box, &QDoubleSpinBox::setValue);

    connect(ui.modulus_of_elasticity_spin_box,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            ui.openGLWidget, &GLWidget::SetModulusOfElasticity);
    connect(ui.openGLWidget, &GLWidget::OnModulusOfElasticityChange,
            ui.modulus_of_elasticity_spin_box, &QDoubleSpinBox::setValue);
}
