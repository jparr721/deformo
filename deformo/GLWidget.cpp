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
    camera_ = std::make_unique<Camera>();
}

GLWidget::~GLWidget() {
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ibo);
    glDeleteBuffers(1, &c_vbo);
    glDeleteVertexArrays(1, &vao);

    delete draw_timer_;
}

void GLWidget::Cleanup() { delete shader_program; }

void GLWidget::Update() {
    input_->Update();

    if (input_->MouseButtonPressed(Qt::RightButton)) {
        camera_->Rotate(-1 * camera_->kRotationSpeed * input_->MouseDelta().x(),
                        camera_->kUp);
        camera_->Rotate(-1 * camera_->kRotationSpeed * input_->MouseDelta().y(),
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

    if (input_->KeyPressed(Qt::Key_Space)) {
        camera_->Reset();
    }

    if (input_->KeyPressed(Qt::Key_G)) {
        simulating_ = !simulating_;
    }

    if (input_->KeyPressed(Qt::Key_M)) {
        render_style = render_style == GL_LINE ? GL_FILL : GL_LINE;
    }
}

void GLWidget::SetCutPlane(float value) {
    mesh->SetCutPlane(value / 100.f);
    emit OnCutPlaneChange(value);
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

    shader_program = new QOpenGLShaderProgram(this);

    BuildMesh();
    BuildPhysicsEngine();

    shader_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
    shader_program->addShaderFromSourceFile(QOpenGLShader::Fragment,
                                            "./core.frag");
    shader_program->link();
    shader_program->bind();

    BuildBuffers();
    shader_program->release();

    // Configure camera matrix positions
    model_loc = shader_program->uniformLocation("m");
    view_loc = shader_program->uniformLocation("v");
    projection_loc = shader_program->uniformLocation("p");

    LogErrors("initializeGL");
}

void GLWidget::paintGL() {
    // Wire
    glPolygonMode(GL_FRONT_AND_BACK, render_style);

    glClear(GL_COLOR_BUFFER_BIT);

    shader_program->bind();

    shader_program->setUniformValue(projection_loc, camera_->Matrix());

    BindVertexAttributeArray(shader_program->programId(), "position", vbo, 3, mesh->positions);
    BindVertexAttributeArray(shader_program->programId(), "color", c_vbo, 3, mesh->colors);

    // Render
    glBindVertexArray(vao);

    glDrawElements(GL_TRIANGLES, mesh->faces_size(), GL_UNSIGNED_INT, nullptr);

    glBindVertexArray(0);

    shader_program->release();

    // Solve at this timestep
    if (simulating_) {
        const auto cycles = utils::stopwatch::time([&] { sim->Solve(); });
        std::cout << "Solver took: " << cycles.count() << " cycles"
                  << std::endl;
    }

    LogErrors("paintGL");
}

void GLWidget::BuildBuffers() {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &c_vbo);

    BindVertexAttributeArray(shader_program->programId(), "position", vbo, 3, mesh->positions);
    BindVertexAttributeArray(shader_program->programId(), "color", c_vbo, 3, mesh->colors);
    BindElementArrayObject(ibo, mesh->faces);
}

void GLWidget::BuildMesh(const float cut_plane) {
    const std::string cdir = std::filesystem::current_path().string();

    const auto ply_path =
        std::filesystem::path(cdir + "/cube.ply");

    mesh = std::make_shared<Mesh>(ply_path.string(), cut_plane);
}

void GLWidget::BuildPhysicsEngine() {
    assert(mesh != nullptr && "MESH IS NOT INITIALIZED");
    std::vector<unsigned int> dynamic_indices;
    utils::FindMaxVertices(dynamic_indices, mesh->positions);

    for (const auto& face_index : dynamic_indices) {
        mesh->colors.segment(face_index, 3) << 0.f, 1.f, 0.f;
    }

    const auto uniform_gravity = Eigen::Vector3f(0.f, -9.81f, 0.f);
    const auto boundary_conditions =
        AssignBoundaryConditionToFixedNodes(dynamic_indices, uniform_gravity);
    sim = std::make_unique<LinearTetrahedral>(210e6, 0.3, 1.f, mesh,
                                              boundary_conditions);
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

void GLWidget::LogErrors(const char* fn) {
    GLenum err;
    for (;;) {
        err = glGetError();

        if (err == GL_NO_ERROR) {
            break;
        }

        std::cerr << "Error in fn: " << fn << ": " << err << std::endl;
    }
}
