#include "Window.h"

Window::Window(QWidget *parent) : QMainWindow(parent) {
  ui.setupUi(this);

  connect(ui.cut_plane_slider, &QSlider::valueChanged, ui.openGLWidget,
          &GLWidget::SetCutPlane);
  connect(ui.openGLWidget, &GLWidget::OnCutPlaneChange, ui.cut_plane_slider,
          &QSlider::setValue);
}
