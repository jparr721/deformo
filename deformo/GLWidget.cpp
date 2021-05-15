#define _USE_MATH_DEFINES

#include "GLWidget.h"

#include <Eigen/Dense>
#include <QVector3D>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

void GLWidget::Cleanup() {}

void GLWidget::initializeGL() {
  initializeOpenGLFunctions();
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);

  // White background
  glClearColor(255.f, 255.f, 255.f, 1.f);

  shader_program = std::make_shared<QOpenGLShaderProgram>(this);

  shader_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
  shader_program->addShaderFromSourceFile(QOpenGLShader::Fragment,
                                          "./core.frag");
  shader_program->link();
  shader_program->bind();

  // Configure camera matrix positions
  model_loc = shader_program->uniformLocation("m");
  view_loc = shader_program->uniformLocation("v");
  projection_loc = shader_program->uniformLocation("p");

  projection.perspective(45.f, 4.f / 3.f, 0.f, 2000.f);
  projection.translate(QVector3D(0.f, 0.f, -5.f));
}

void GLWidget::paintGL() {
  // Solve at this timestep
  // sim->Solve();

  // Integrate vars to next step
  // sim->Integrate();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  shader_program->bind();

  shader_program->setUniformValue(projection_loc, projection);

  mesh->Render();

  shader_program->release();
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
