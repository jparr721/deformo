#define _USE_MATH_DEFINES

#include "GLWidget.h"

#include <Eigen/Dense>
#include <QVector3D>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

static std::vector<Vertex> triangle = {
    Vertex(QVector3D(0.00f, 0.75f, 0.f), QVector3D(1.0f, 0.0f, 0.f)),
    Vertex(QVector3D(-0.75f, -0.75f, 0.f), QVector3D(0.0f, 1.0f, 0.f)),
    Vertex(QVector3D(0.75f, -0.75f, 0.f), QVector3D(0.0f, 0.0f, 0.f))};

void GLWidget::Cleanup() { delete program; }

void GLWidget::initializeGL() {
  initializeOpenGLFunctions();
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);

  // White background
  glClearColor(255.f, 255.f, 255.f, 1.f);

  Eigen::VectorXd displacements(12);
   displacements << 0, 0, 0.5, 0, 0.5, 0.25, 0, 0, 0.5, 0.25, 0, 0.25;

  mesh = std::make_shared<Mesh>(displacements, true);
  sim = std::make_unique<Simulation>(1., 210e6, 0.3, mesh,
                                     std::vector<BoundaryCondition>{
                                         {
                                             2,
                                             xyz{9.375, 0},
                                         },
                                         {
                                             3,
                                             xy{9.375, 0},
                                         },
                                     });

  program = new QOpenGLShaderProgram(this);

  program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
  program->addShaderFromSourceFile(QOpenGLShader::Fragment, "./core.frag");
  program->link();
  program->bind();

  // Configure camera matrix positions
  model_loc = program->uniformLocation("m");
  view_loc = program->uniformLocation("v");
  projection_loc = program->uniformLocation("p");

  projection.perspective(45.f, 4.f / 3.f, 0.f, 2000.f);
  projection.translate(QVector3D(0.f, 0.f, -5.f));

  // Passthrough for data to GPU
  mesh->Initialize(vbo);

  // Create vao
  vao.create();
  vao.bind();
  program->enableAttributeArray(0);  // Vertices
  program->enableAttributeArray(1);  // Colors
  program->setAttributeBuffer(0, GL_FLOAT, Vertex::PositionOffset(),
                              Vertex::PositionSize(), Vertex::Stride());
  program->setAttributeBuffer(1, GL_FLOAT, Vertex::ColorOffset(),
                              Vertex::ColorSize(), Vertex::Stride());

  vao.release();
  vbo.release();
  program->release();
}

void GLWidget::paintGL() {
  // Solve at this timestep
  // sim->Solve();

  // Integrate vars to next step
  // sim->Integrate();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  program->bind();

  program->setUniformValue(projection_loc, projection);

  mesh->Render(vbo, vao);

  program->release();
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
