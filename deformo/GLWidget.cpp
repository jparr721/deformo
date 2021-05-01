#include "GLWidget.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "Mesh.h"
#include "Simulation.h"

void GLWidget::Cleanup() { delete program_id; }

void GLWidget::initializeGL() {
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);

  // White background
  glClearColor(255.f, 255.f, 255.f, 1.f);

  Eigen::VectorXd displacements(12);
  std::vector<Eigen::Vector3f> triangle{
      {Eigen::Vector3f(1.f, -1.f, 0.f)},  {Eigen::Vector3f(1.f, 0.f, 0.f)},
      {Eigen::Vector3f(0.f, 1.f, 0.f)},   {Eigen::Vector3f(0.f, 1.f, 0.f)},
      {Eigen::Vector3f(-1.f, -1.f, 0.f)}, {Eigen::Vector3f(0.f, 0.f, 1.f)},
  };

  displacements << 0, 0, 0.5, 0, 0.5, 0.25, 0, 0, 0.5, 0.25, 0, 0.25;

  const auto mesh = std::make_shared<Mesh<2>>(displacements);
  const auto sim = std::make_unique<Simulation>(1., 210e6, 0.3, mesh,
                                                std::vector<BoundaryCondition>{
                                                    {
                                                        2,
                                                        xy{9.375, 0},
                                                    },
                                                    {
                                                        3,
                                                        xy{9.375, 0},
                                                    },
                                                });
  sim->Solve();

  initializeOpenGLFunctions();

  program_id = new QOpenGLShaderProgram(this);

  program_id->addShaderFromSourceFile(QOpenGLShader::Vertex, "./core.vs");
  program_id->addShaderFromSourceFile(QOpenGLShader::Fragment, "./core.frag");
  program_id->link();
  program_id->bind();

  // Passthrough for data to GPU
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(triangle.data(), triangle.size());

  // Create vao
  vao.create();
  vao.bind();
  program_id->enableAttributeArray(0);  // Vertices
  program_id->enableAttributeArray(1);  // Colors
  program_id->setAttributeBuffer(0, GL_FLOAT, Vertex<2>::PositionOffset(),
                                 Vertex<2>::Size(), Vertex<2>::Stride());
  program_id->setAttributeBuffer(1, GL_FLOAT, Vertex<2>::ColorOffset(),
                                 Vertex<2>::Size(), Vertex<2>::Stride());

  vbo.release();
  vao.release();
  program_id->release();
}

void GLWidget::paintGL() {
  fov.perspective(45.f, 4.f / 3.f, 0.1f, 2000.f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  program_id->bind();

  program_id->setUniformValue(matrix_uniform, fov);

  glVertexAttribPointer(position, 3, GL_FLOAT, GL_FALSE, 0,
                        static_cast<void*>(triangle.data()));
  glVertexAttribPointer(color, 3, GL_FLOAT, GL_FALSE, 0,
                        static_cast<void*>(colors.data()));

  glEnableVertexAttribArray(position);
  glEnableVertexAttribArray(color);

  glDrawArrays(GL_TRIANGLES, 0, triangle.size());

  glDisableVertexAttribArray(position);
  glDisableVertexAttribArray(color);

  program_id->release();
}
