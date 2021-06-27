#include "Window.h"

Window::Window(QWidget* parent) : QMainWindow(parent) {
    ui.setupUi(this);
    controller_ = std::make_shared<WindowController>(ui, "cube.ply");
    ui.openGLWidget->SetController(controller_);
}
