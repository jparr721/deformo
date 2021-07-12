#include "Window.h"

Window::Window(QWidget* parent) : QMainWindow(parent) {
    ui.setupUi(this);
    controller_ = std::make_shared<WindowController>(ui);
    ui.openGLWidget->SetController(controller_);
    ui.openGLWidget_2->SetController(controller_);
}
