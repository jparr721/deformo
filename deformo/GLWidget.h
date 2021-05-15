#pragma once

#include <Eigen/Dense>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <memory>
#include <string>

#include "Mesh.h"
#include "Simulation.h"

void Perspective(Eigen::Matrix4d& camera, float vertical_angle,
                 float aspect_ratio, float near_plane, float far_plane);

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
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

  GLWidget() {}
  GLWidget(QWidget* parent = nullptr) : QOpenGLWidget(parent) {}

  void resizeGL(int width, int height) override;

 public slots:
  void Cleanup();

 protected:
  void initializeGL() override;
  void paintGL() override;
};
