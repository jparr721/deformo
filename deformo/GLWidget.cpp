#define _USE_MATH_DEFINES

#include "GLWidget.h"

#include <Eigen/Dense>
#include <QVector3D>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "Mesh.h"
#include "Simulation.h"

static VVertex tt[] = {
    VVertex(QVector3D(0.00f, 0.75f, 1.0f), QVector3D(1.0f, 0.0f, 0.0f)),
    VVertex(QVector3D(-0.75f, -0.75f, 1.0f), QVector3D(0.0f, 1.0f, 0.0f)),
    VVertex(QVector3D(0.75f, -0.75f, 1.0f), QVector3D(0.0f, 0.0f, 1.0f))};

static std::vector<Vertex> triangle = {
    Vertex(QVector3D(0.00f, 0.75f, 1.0f), QVector3D(1.0f, 0.0f, 0.0f)),
    Vertex(QVector3D(-0.75f, -0.75f, 1.0f), QVector3D(0.0f, 1.0f, 0.0f)),
    Vertex(QVector3D(0.75f, -0.75f, 1.0f), QVector3D(0.0f, 0.0f, 1.0f))};

void GLWidget::Cleanup() { delete program; }

void GLWidget::initializeGL() {
  initializeOpenGLFunctions();
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);

  std::cout << tt[0].positionOffset() << std::endl;
  std::cout << tt[0].colorOffset() << std::endl;
  std::cout << tt[0].stride() << std::endl;
  std::cout << sizeof(tt) / sizeof(tt[0]) << std::endl;

  std::cout << std::endl;
  std::cout << triangle[0].PositionOffset() << std::endl;
  std::cout << triangle[0].ColorOffset() << std::endl;
  std::cout << triangle[0].Stride() << std::endl;
  std::cout << triangle.size() * sizeof(Vertex) << std::endl;

  // White background
  // glClearColor(255.f, 255.f, 255.f, 1.f);
  glClearColor(0.f, 0.f, 0.f, 1.f);

  // Eigen::VectorXd displacements(12);

  // displacements << 0, 0, 0.5, 0, 0.5, 0.25, 0, 0, 0.5, 0.25, 0, 0.25;

  // const auto mesh = std::make_shared<Mesh>(displacements);
  // const auto sim = std::make_unique<Simulation>(1., 210e6, 0.3, mesh,
  //                                              std::vector<BoundaryCondition>{
  //                                                  {
  //                                                      2,
  //                                                      xy{9.375, 0},
  //                                                  },
  //                                                  {
  //                                                      3,
  //                                                      xy{9.375, 0},
  //                                                  },
  //                                              });
  // sim->Solve();

  program = new QOpenGLShaderProgram(this);

  program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
  program->addShaderFromSourceFile(QOpenGLShader::Fragment, "./core.frag");
  program->link();
  program->bind();

  // Configure camera matrix positions
  model_loc = program->uniformLocation("m");
  view_loc = program->uniformLocation("v");
  projection_loc = program->uniformLocation("p");

  Perspective(projection, 45.f, 4.f / 3.f, 0.1f, 2000.f);

  // Passthrough for data to GPU
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(triangle.data(), triangle.size() * sizeof(Vertex));

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
  //QMatrix4x4 proj;
  //proj.perspective(45.f, 4.f / 3.f, 0.1f, 2000.f);
  //proj.translate(0, 0, 5);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  program->bind();

  //program->setUniformValue(projection_loc, proj);
  // triangle[0]
  //    .position = (QVector3D(triangle[0].position.x() + 0.1,
  //    triangle[0].position.y(), triangle[0].position.z()));

  // vbo.bind();
  // vbo.write(Vertex::PositionOffset(), static_cast<void*>(triangle.data()),
  //          triangle.size() * sizeof(QVector3D));
  // vbo.release();

  vao.bind();
  glDrawArrays(GL_TRIANGLES, 0, triangle.size());
  vao.release();

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
