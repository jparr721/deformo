#pragma once

#include <QMatrix4x4>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <string>

std::string ReadFileToString(const std::string& path);

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
  QOpenGLShaderProgram* program_id;
  QMatrix4x4 fov;

  GLint position = 0;
  GLint color = 0;
  GLint matrix_uniform = 0;

  GLWidget(QWidget* parent = nullptr) : QOpenGLWidget(parent) {}

 public slots:
  void Cleanup();

 protected:
  void initializeGL() override;
  void paintGL() override;
};
