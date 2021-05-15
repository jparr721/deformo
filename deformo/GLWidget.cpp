#define _USE_MATH_DEFINES

#include "GLWidget.h"

#include <Eigen/Dense>
#include <QVector3D>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Loader.h"

GLWidget::GLWidget(QWidget* parent) : QOpenGLWidget(parent) {
  setFocusPolicy(Qt::ClickFocus);
  input = std::make_unique<Input>();
  camera = std::make_unique<Camera>();
}
void GLWidget::Cleanup() { delete shader_program; }

void GLWidget::Update() {
  if (input->KeyPressed(Qt::Key_W)) {
    projection.translate(camera->kForward);
  }

  if (input->KeyPressed(Qt::Key_S)) {
    projection.translate(-camera->kForward);
  }
}

void GLWidget::initializeGL() {
  initializeOpenGLFunctions();
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);
  connect(this, &QOpenGLWidget::frameSwapped, this, &GLWidget::Update);

  // White background
  glClearColor(255.f, 255.f, 255.f, 1.f);

  shader_program = new QOpenGLShaderProgram(this);

  // Initialize Mesh ====================
  Eigen::MatrixXf V;
  Eigen::MatrixXf F;
  Eigen::MatrixXf T;

  const std::string cdir = std::filesystem::current_path().string();

  const std::filesystem::path node_path =
      std::filesystem::path(cdir + "/square.1.node");
  const std::filesystem::path ele_path =
      std::filesystem::path(cdir + "/square.1.ele");
  const std::filesystem::path face_path =
      std::filesystem::path(cdir + "/square.1.face");

  loader::ReadTetgenVertexFile(V, node_path.string());
  loader::ReadTetgenFaceFile(F, face_path.string());
  loader::ReadTetgenEleFile(T, ele_path.string());

  mesh = std::make_shared<Mesh>(V, F, T);

  // ====================
  shader_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
  shader_program->addShaderFromSourceFile(QOpenGLShader::Fragment,
                                          "./core.frag");
  shader_program->link();
  shader_program->bind();

  vbo = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(mesh->data(), mesh->size_bytes());

  //ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
  //ibo.create();
  //ibo.bind();
  //ibo.allocate(mesh->faces.data(), mesh->faces.size() * sizeof(float));
  //// Indices won't be changing
  //ibo.setUsagePattern(QOpenGLBuffer::StaticDraw);

  vao.create();
  vao.bind();

  // Vertices
  shader_program->enableAttributeArray(0);

  // Colors
  shader_program->enableAttributeArray(1);

  shader_program->setAttributeBuffer(0, GL_FLOAT, 0, 3, sizeof(Vertex));
  shader_program->setAttributeBuffer(1, GL_FLOAT, offsetof(Vertex, color), 3, sizeof(Vertex));
  //shader_program->setAttributeArray(0, GL_FLOAT, mesh->vertices.data(), 3);
  //shader_program->setAttributeArray(1, GL_FLOAT, mesh->colors.data(), 3);

  vao.release();
  vbo.release();
  //ibo.release();
  shader_program->release();

  // Configure camera matrix positions
  model_loc = shader_program->uniformLocation("m");
  view_loc = shader_program->uniformLocation("v");
  projection_loc = shader_program->uniformLocation("p");

  projection.perspective(45.f, 4.f / 3.f, 1.f, 2000.f);
  projection.translate(QVector3D(0.f, 0.f, -5.f));
}

void GLWidget::paintGL() {
  // Solve at this timestep
  // sim->Solve();

  // Integrate vars to next step
  // sim->Integrate();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  shader_program->bind();

  //shader_program->setUniformValue(view_loc, camera->Matrix());
  shader_program->setUniformValue(projection_loc, projection);

  // Add updated vertex coordinates
  vbo.bind();
  vbo.write(0, mesh->data(), mesh->size_bytes());
  vbo.release();

  // Render
  vao.bind();
  //glDrawElements(GL_TRIANGLES, mesh->size(), GL_FLOAT, 0);
  glDrawArrays(GL_TRIANGLES, 0, mesh->size());
  vao.release();

  shader_program->release();
}

void GLWidget::keyPressEvent(QKeyEvent* event) {
  input->RegisterKeyPress(event->key());
}

void GLWidget::keyReleaseEvent(QKeyEvent* event) {
  input->RegisterKeyRelease(event->key());
}

void GLWidget::resizeGL(int width, int height) {}

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
