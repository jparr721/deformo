#pragma once

#include <QMatrix4x4>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <string>

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
	 // Global state vars
  QOpenGLBuffer vbo;
  QOpenGLVertexArrayObject vao;
  QOpenGLShaderProgram* program_id;

  // Camera vars
  int model;
  int world;
  int 
  QMatrix4x4 projection;

  GLWidget(QWidget* parent = nullptr) : QOpenGLWidget(parent) {}

 public slots:
  void Cleanup();

 protected:
  void initializeGL() override;
  void paintGL() override;
};
