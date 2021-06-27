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
    camera_ = std::make_shared<Camera>();
}

GLWidget::~GLWidget() {
    delete draw_timer_;
}

void GLWidget::Update() {
    input_->Update();

    if (input_->MouseButtonPressed(Qt::RightButton)) {
        camera_->Rotate(camera_->kRotationSpeed * input_->MouseDelta().x(),
                        camera_->kUp);
        camera_->Rotate(camera_->kRotationSpeed * input_->MouseDelta().y(),
                        camera_->rotation.rotatedVector(camera_->kRight));
    }

    if (input_->KeyPressed(Qt::Key_W)) {
        camera_->Forward();
    }

    if (input_->KeyPressed(Qt::Key_S)) {
        camera_->Backward();
    }

    if (input_->KeyPressed(Qt::Key_A)) {
        camera_->Left();
    }

    if (input_->KeyPressed(Qt::Key_D)) {
        camera_->Right();
    }

    if (input_->KeyPressed(Qt::Key_Q)) {
        camera_->Up();
    }

    if (input_->KeyPressed(Qt::Key_E)) {
        camera_->Down();
    }

    if (input_->KeyPressed(Qt::Key_R)) {
        controller_->Reset();
    }

    if (input_->KeyPressed(Qt::Key_Space)) {
        camera_->Reset();
    }

    if (input_->KeyPressed(Qt::Key_G)) {
        simulating_ = !simulating_;
    }

    if (input_->KeyPressed(Qt::Key_M)) {
        renderer->SetRenderMode();
    }
}

void GLWidget::initializeGL() {
    initializeOpenGLFunctions();
    if (const auto code = glewInit(); code != GLEW_OK) {
        std::cerr << glewGetErrorString(code) << std::endl;
        exit(EXIT_FAILURE);
    }
    connect(this, &QOpenGLWidget::frameSwapped, this, &GLWidget::Update);

    // Face Culling
    glEnable(GL_DEPTH);

    // White background
    glClearColor(255.f, 255.f, 255.f, 1.f);

    shader_program = std::make_shared<ShaderProgram>();

    renderer = std::make_unique<Renderer>(controller_->mesh, shader_program, camera_);

    LogErrors("initializeGL");
}

void GLWidget::paintGL() {
    renderer->Render();

    if (simulating_) {
        controller_->StepForward();
    }

    LogErrors("paintGL");
}

void GLWidget::keyPressEvent(QKeyEvent* event) {
    input_->RegisterKeyPress(event->key());
}

void GLWidget::keyReleaseEvent(QKeyEvent* event) {
    input_->RegisterKeyRelease(event->key());
}

void GLWidget::mousePressEvent(QMouseEvent* event) {
    input_->RegisterMouseButtonPress(event->button());
}

void GLWidget::mouseReleaseEvent(QMouseEvent* event) {
    input_->RegisterMouseButtonRelease(event->button());
}

void GLWidget::resizeGL(int width, int height) {}
