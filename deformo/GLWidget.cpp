#define _USE_MATH_DEFINES

#include "GLWidget.h"
#include "Utils.h"

#include <Eigen/Dense>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

GLWidget::GLWidget(QWidget* parent) : QOpenGLWidget(parent) {
    setFocusPolicy(Qt::ClickFocus);

    draw_timer = new QTimer(this);
    connect(draw_timer, &QTimer::timeout, this,
            QOverload<>::of(&GLWidget::update));
    connect(draw_timer, &QTimer::timeout, this,
            QOverload<>::of(&GLWidget::Update));
    draw_timer->start(30);

    input = std::make_unique<Input>();
    camera = std::make_unique<Camera>();
}

GLWidget::~GLWidget() {
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ibo);
    glDeleteBuffers(1, &c_vbo);
    vao.destroy();

    delete draw_timer;
}

void GLWidget::Cleanup() { delete shader_program; }

void GLWidget::Update() {
    input->Update();

    if (input->MouseButtonPressed(Qt::RightButton)) {
        camera->Rotate(-1 * camera->kRotationSpeed * input->MouseDelta().x(),
                       camera->kUp);
        camera->Rotate(-1 * camera->kRotationSpeed * input->MouseDelta().y(),
                       camera->rotation.rotatedVector(camera->kRight));
    }

    if (input->KeyPressed(Qt::Key_W)) {
        camera->Forward();
    }

    if (input->KeyPressed(Qt::Key_S)) {
        camera->Backward();
    }

    if (input->KeyPressed(Qt::Key_A)) {
        camera->Left();
    }

    if (input->KeyPressed(Qt::Key_D)) {
        camera->Right();
    }

    if (input->KeyPressed(Qt::Key_Q)) {
        camera->Up();
    }

    if (input->KeyPressed(Qt::Key_E)) {
        camera->Down();
    }

    if (input->KeyPressed(Qt::Key_Space)) {
        camera->Reset();
    }

    if (input->KeyPressed(Qt::Key_M)) {
        if (render_style == GL_LINE) {
            render_style = GL_FILL;
        } else {
            render_style = GL_LINE;
        }
    }
}

void GLWidget::SetCutPlane(float value) {
    mesh->SetCutPlane(value / 100.f);
    emit OnCutPlaneChange(value);
}

void GLWidget::initializeGL() {
    initializeOpenGLFunctions();
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

    // Solve at this timestep
    while (simulating) {
        const auto cycles = utils::stopwatch::time([&] { sim->Solve(); });
        std::cout << "Solver took: " << cycles.count() << " cycles"
                  << std::endl;
    }

    glClear(GL_COLOR_BUFFER_BIT);

    shader_program->bind();

    shader_program->setUniformValue(projection_loc, camera->Matrix());

    // Add updated vertex coordinates
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, mesh->size_bytes(), mesh->data(),
                 GL_DYNAMIC_DRAW);

    // Add updated colors
    glBindBuffer(GL_ARRAY_BUFFER, c_vbo);
    glBufferData(GL_ARRAY_BUFFER, mesh->colors_size_bytes(),
                 mesh->colors_data(), GL_DYNAMIC_DRAW);

    // Render
    vao.bind();

    glDrawElements(GL_TRIANGLES, mesh->faces_size(), GL_UNSIGNED_INT, 0);

    vao.release();

    shader_program->release();
    LogErrors("paintGL");
}

void GLWidget::BuildBuffers() {
    vao.create();
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &c_vbo);
    glGenBuffers(1, &ibo);

    vao.bind();

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, mesh->size_bytes(), mesh->data(),
                 GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);

    glBindBuffer(GL_ARRAY_BUFFER, c_vbo);
    glBufferData(GL_ARRAY_BUFFER, mesh->colors_size_bytes(),
                 mesh->colors_data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->faces_size_bytes(),
                 mesh->faces_data(), GL_DYNAMIC_DRAW);
}

void GLWidget::BuildMesh(const float cut_plane) {
    const std::string cdir = std::filesystem::current_path().string();

    const std::filesystem::path ply_path =
        std::filesystem::path(cdir + "/cube.ply");

    mesh = std::make_shared<Mesh>(ply_path.string(), cut_plane);
}

void GLWidget::BuildPhysicsEngine() {
    assert(mesh != nullptr && "MESH IS NOT INITIALIZED");
    std::vector<unsigned int> dynamic_indices;
    utils::FindMaxVertices(dynamic_indices, mesh->positions);
    const Eigen::Vector3f uniform_gravity = Eigen::Vector3f(0.f, 9.81f, 0.f);
    const auto boundary_conditions =
        AssignBoundaryConditionToFixedNodes(dynamic_indices, uniform_gravity);
    sim = std::make_unique<LinearTetrahedral>(210e6, 0.3, 1.f, mesh,
                                              boundary_conditions);
}

void GLWidget::keyPressEvent(QKeyEvent* event) {
    input->RegisterKeyPress(event->key());
}

void GLWidget::keyReleaseEvent(QKeyEvent* event) {
    input->RegisterKeyRelease(event->key());
}

void GLWidget::mousePressEvent(QMouseEvent* event) {
    input->RegisterMouseButtonPress(event->button());
}

void GLWidget::mouseReleaseEvent(QMouseEvent* event) {
    input->RegisterMouseButtonRelease(event->button());
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

void Perspective(Eigen::Matrix4d& camera, float vertical_angle,
                 float aspect_ratio, float near_plane, float far_plane) {
    assert(far_plane > near_plane &&
           "NEAR PLANE CANNOT BE GREATER THAN OR EQUAL TO FAR PLANE");

    Eigen::Matrix4d m;
    const float vertical_angle_rad = vertical_angle * (M_PI / 180);
    const float sine = std::sin(vertical_angle_rad);

    if (sine == 0.f) {
        return;
    }

    const float cotan = std::cos(vertical_angle_rad) / sine;
    const float clip = far_plane - near_plane;

    m(0, 0) = cotan / aspect_ratio;
    m(1, 1) = cotan;
    m(2, 2) = -(near_plane + far_plane) / clip;
    m(2, 3) = -1.f;
    m(3, 2) = -(2.0f * near_plane * far_plane) / clip;

    camera *= m;
}
