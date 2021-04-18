#include "GLWidget.h"

#include <Eigen/Dense>
#include <fstream>
#include <sstream>

void GLWidget::Cleanup() { delete program_id; }

void GLWidget::initializeGL() {
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this,
          &GLWidget::Cleanup);
  initializeOpenGLFunctions();

  program_id = new QOpenGLShaderProgram(this);

  const auto vertex_shader = ReadFileToString("./core.vs");
  const auto fragment_shader = ReadFileToString("./core.frag");

  program_id->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                      vertex_shader.data());
  program_id->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                      fragment_shader.data());
  program_id->link();

  position = program_id->attributeLocation("position_attribute");
  Q_ASSERT(position > 0);
  color = program_id->attributeLocation("color_attribute");
  Q_ASSERT(color > 0);
  matrix_uniform = program_id->uniformLocation("projection_matrix");
  Q_ASSERT(matrix_uniform > 0);
}

void GLWidget::paintGL() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  program_id->bind();

  fov.perspective(45.f, 4.f / 3.f, 0.1f, 2000.f);
  fov.translate(0, 0, -5);

  program_id->setUniformValue(matrix_uniform, fov);

  std::vector<Eigen::Vector3f> triangle{
      {Eigen::Vector3f(1.f, -1.f, 0.f)},
      {Eigen::Vector3f(0.f, 1.f, 0.f)},
      {Eigen::Vector3f(-1.f, -1.f, 0.f)},
  };
  std::vector<Eigen::Vector3f> colors{
      {Eigen::Vector3f(1.f, 0.f, 0.f)},
      {Eigen::Vector3f(0.f, 1.f, 0.f)},
      {Eigen::Vector3f(0.f, 0.f, 1.f)},
  };

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

std::string ReadFileToString(const std::string& path) {
  std::ostringstream sstr;
  const auto stream = std::ifstream{path};
  sstr << stream.rdbuf();
  return sstr.str();
}
