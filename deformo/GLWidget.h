#pragma once

#include <Eigen/Dense>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QKeyEvent>
#include <memory>
#include <string>

#include "Camera.h"
#include "Input.h"
#include "Mesh.h"
#include "Simulation.h"

void Perspective(Eigen::Matrix4d& camera, float vertical_angle,
                 float aspect_ratio, float near_plane, float far_plane);

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
  // Vertex Buffer
  QOpenGLBuffer vbo;

  // Index Buffer
  QOpenGLBuffer ibo;

  // Vertex Array Object
  QOpenGLVertexArrayObject vao;

  QOpenGLShaderProgram* shader_program;
  // Camera Shader Locations
  int model_loc = 0;
  int view_loc = 0;
  int projection_loc = 0;

  // Camera Shader Values
  QMatrix4x4 model;
  QMatrix4x4 view;
  QMatrix4x4 projection;

  // Mesh object for 3D geometry
  std::shared_ptr<Mesh> mesh;

  // Simulation Object
  std::unique_ptr<Simulation> sim;

  GLWidget(QWidget* parent = nullptr);

  void resizeGL(int width, int height) override;

 public slots:
  void Cleanup();
  void Update();

 protected:
  void initializeGL() override;
  void paintGL() override;

  // Keyboard Shenanigans
  void keyPressEvent(QKeyEvent* event) override;
  void keyReleaseEvent(QKeyEvent* event) override;

 private:
  std::unique_ptr<Input> input;
  std::unique_ptr<Camera> camera;
};
