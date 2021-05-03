#pragma once

#include <Eigen/Dense>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <string>
#include "Vertex.h"

void Perspective(Eigen::Matrix4d& camera, float vertical_angle,
                 float aspect_ratio, float near_plane, float far_plane);

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
  // const std::vector<Eigen::Vector3f> triangle{
  //    {Eigen::Vector3f(1.f, -1.f, 0.f)},  {Eigen::Vector3f(1.f, 0.f, 0.f)},
  //    {Eigen::Vector3f(0.f, 1.f, 0.f)},   {Eigen::Vector3f(0.f, 1.f, 0.f)},
  //    {Eigen::Vector3f(-1.f, -1.f, 0.f)}, {Eigen::Vector3f(0.f, 0.f, 1.f)},
  //};

  //const float triangle[18] = {1.f, -1.f, 0.f, 1.f,  0.f,  0.f, 0.f, 1.f, 0.f,
  //                          0.f, 1.f,  0.f, -1.f, -1.f, 0.f, 0.f, 0.f, 1.f};

  // Global state vars
  QOpenGLBuffer vbo;
  QOpenGLVertexArrayObject vao;
  QOpenGLShaderProgram* program;

  // Camera Shader Locations
  int model_loc;
  int view_loc;
  int projection_loc;

  // Camera Shader Values
  Eigen::Matrix4d model;
  Eigen::Matrix4d view;
  Eigen::Matrix4d projection;

  GLWidget(QWidget* parent = nullptr) : QOpenGLWidget(parent) {}

  void resizeGL(int width, int height) override;

 public slots:
  void Cleanup();

 protected:
  void initializeGL() override;
  void paintGL() override;
};
