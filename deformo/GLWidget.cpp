#define _USE_MATH_DEFINES

#include "GLWidget.h"
#include "Utils.h"

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
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ibo);
    glDeleteBuffers(1, &c_vbo);
    glDeleteVertexArrays(1, &vao);

    delete draw_timer_;
}

void GLWidget::Cleanup() {}

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
        Reset();
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
    connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
            &GLWidget::Cleanup);
    connect(this, &QOpenGLWidget::frameSwapped, this, &GLWidget::Update);

    // Face Culling
    glEnable(GL_DEPTH);

    // White background
    glClearColor(255.f, 255.f, 255.f, 1.f);

    shader_program = std::make_shared<ShaderProgram>();

    BuildMesh();
    BuildPhysicsEngine();

    renderer = std::make_unique<Renderer>(mesh, shader_program, camera_);

    LogErrors("initializeGL");
}

void GLWidget::paintGL() {
    renderer->Render();

    if (simulating_) {
        sim->Solve();
    }

    LogErrors("paintGL");
}

void GLWidget::BuildMesh(const float cut_plane) {
    const std::string cdir = std::filesystem::current_path().string();

    const auto ply_path = std::filesystem::path(cdir + "/cube.ply");

    mesh = std::make_shared<Mesh>(ply_path.string(), cut_plane);
}

void GLWidget::BuildPhysicsEngine() {
    assert(mesh != nullptr && "MESH IS NOT INITIALIZED");

    const auto uniform_gravity = Eigen::Vector3f(0.f, -9.81f, 0.f);
    std::vector<unsigned int> dynamic_indices;
    utils::FindMaxVertices(dynamic_indices, mesh->positions);

    for (const auto& face_index : dynamic_indices) {
        mesh->colors.segment((face_index / 3) * 4, 4)
            << Mesh::kMeshDefaultSelectedColor;
    }
    const auto boundary_conditions =
        AssignBoundaryConditionToFixedNodes(dynamic_indices, uniform_gravity);
    sim = std::make_unique<LinearTetrahedral>(1000, 0.3, 1.f, mesh,
                                              boundary_conditions);
}

void GLWidget::Reset() {
    BuildMesh();
    BuildPhysicsEngine();
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

// SLOTS ====================
void GLWidget::SetSliceAxis(const QString& value) {
    emit OnSliceAxisChange(value);
}

void GLWidget::SetSliceValue(float value) {
    mesh->SetCutPlane(value / 100.f);
    emit OnSliceValueChange(value);
}

void GLWidget::SetNodalMass(float value) { emit OnNodalMassChange(value); }

void GLWidget::SetPoissonsRatio(double value) {
    emit OnPoissonsRatioChange(value);
}

void GLWidget::SetYoungsModulus(double value) {
    emit OnYoungsModulusChange(value);
}

void GLWidget::SetTimestepSize(double value) {
    emit OnTimestepSizeChange(value);
}

void GLWidget::SetRayleighLambda(double value) {
    emit OnRayleighLambdaChange(value);
}

void GLWidget::SetRayleighMu(double value) { emit OnRayleighMuChange(value); }

void GLWidget::RunSimulationButtonPressed() {
    std::cout << "Pressed" << std::endl;
}

void GLWidget::SetTetgenFlags(const QString& value) {
    emit OnTetgenFlagsChange(value);
}
