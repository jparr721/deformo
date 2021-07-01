#define _USE_MATH_DEFINES

#include "GLWidget.h"
#include "QTUtils.h"
#include "Utils.h"
#include "WindowController.h"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>

GLWidget::GLWidget(QWidget* parent) : QOpenGLWidget(parent) {
    setFocusPolicy(Qt::ClickFocus);

    draw_timer_ = new QTimer(this);
    connect(draw_timer_, &QTimer::timeout, this,
            QOverload<>::of(&GLWidget::update));
    connect(draw_timer_, &QTimer::timeout, this,
            QOverload<>::of(&GLWidget::Update));
    draw_timer_->start(30);

    input_ = std::make_unique<Input>();
    camera_ = std::make_shared<Camera<Real>>();
}

GLWidget::~GLWidget() { delete draw_timer_; }

void GLWidget::Update() {}

void GLWidget::initializeGL() {
    initializeOpenGLFunctions();
    if (const auto code = glewInit(); code != GLEW_OK) {
        std::cerr << glewGetErrorString(code) << std::endl;
        exit(EXIT_FAILURE);
    }
    connect(this, &QOpenGLWidget::frameSwapped, this, &GLWidget::Update);

    glEnable(GL_DEPTH_TEST);

    // White background
    glClearColor(255.f, 255.f, 255.f, 1.f);

    shader_program = std::make_shared<ShaderProgram>();

    renderer =
        std::make_shared<Renderer>(controller_->mesh, shader_program, camera_);
    controller_->SetRenderer(renderer);

    LogErrors("initializeGL");
}

void GLWidget::paintGL() {
    renderer->Render();

    if (controller_->IsSimulating()) {
        controller_->StepForward();
    }

    LogErrors("paintGL");
}

void GLWidget::mousePressEvent(QMouseEvent* event) {
    input_->MousePressEvent(event, camera_);
}

void GLWidget::mouseReleaseEvent(QMouseEvent* event) {
    input_->MouseReleaseEvent(event, camera_);
}

void GLWidget::mouseMoveEvent(QMouseEvent* event) {
    input_->MouseMoveEvent(event, camera_);
}

void GLWidget::wheelEvent(QWheelEvent* event) {
    input_->WheelEvent(event, camera_);
}

void GLWidget::resizeGL(int width, int height) {
    renderer->Resize(width, height);
}
