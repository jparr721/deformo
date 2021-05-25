#pragma once

#include <Eigen/Dense>
#include <QKeyEvent>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QTimer>
#include <memory>
#include <string>

#include "Camera.h"
#include "Input.h"
#include "Mesh.h"
#include "LinearTetrahedral.h"

void Perspective(Eigen::Matrix4d& camera, float vertical_angle,
                 float aspect_ratio, float near_plane, float far_plane);

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
  Q_OBJECT

 public:
  unsigned int vbo;
  unsigned int c_vbo;
  unsigned int ibo;

  // Vertex Array Object
  QOpenGLVertexArrayObject vao;

  QOpenGLShaderProgram* shader_program;

  // Toggleable Wire Mesh
  GLenum render_style = GL_LINE;

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

  // LinearTetrahedral Object
  std::unique_ptr<LinearTetrahedral> sim;

  GLWidget(QWidget* parent = nullptr);
  ~GLWidget();

  void resizeGL(int width, int height) override;

 public slots:
  void Cleanup();
  void Update();

  void SetCutPlane(float value);

 signals:
  void OnCutPlaneChange(float value);

 protected:
  void initializeGL() override;
  void paintGL() override;

  // Keyboard Shenanigans
  void keyPressEvent(QKeyEvent* event) override;
  void keyReleaseEvent(QKeyEvent* event) override;

  // Mouse Shenanigans
  void mousePressEvent(QMouseEvent* event) override;
  void mouseReleaseEvent(QMouseEvent* event) override;

 private:
  QTimer* draw_timer;
  std::unique_ptr<Input> input;
  std::unique_ptr<Camera> camera;

  void LogErrors(const char* fn);
  void BuildBuffers();
  void BuildMesh(float cut_plane = Mesh::kNoCutPlane);
};
