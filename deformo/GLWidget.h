#pragma once

#include "Vbo.h"
#include "Vao.h"

#include <Eigen/Dense>
#include <QKeyEvent>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <QTimer>
#include <memory>

#include "Camera.h"
#include "Input.h"
#include "LinearTetrahedral.h"
#include "Mesh.h"

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

  public:
    unsigned int vao = 0;
    unsigned int vbo = 0;
    unsigned int c_vbo = 0;
    unsigned int ibo = 0;

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

    explicit GLWidget(QWidget* parent = nullptr);
    ~GLWidget() override;

    void resizeGL(int width, int height) override;

  public slots:
    void Cleanup();
    void Update();

    void SetCutPlane(float value);
    void SetInteractiveModeToggle(int state);

    void SetPoissonsRatio(double value);
    void SetModulusOfElasticity(double value);

    void SetTetgenFlags(const QString& value);

  signals:
    void OnCutPlaneChange(float value);
    void OnInteractiveModeToggled(int state);

    void OnPoissonsRatioChange(double value);
    void OnModulusOfElasticityChange(double value);

    void OnTetgenFlagsChange(const QString& value);

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
    QTimer* draw_timer_;
    std::unique_ptr<Input> input_;
    std::unique_ptr<Camera> camera_;
    bool simulating_ = false;

    void LogErrors(const char* fn);
    void BuildBuffers();
    void BuildMesh(float cut_plane = Mesh::kNoCutPlane);
    void BuildPhysicsEngine();

    void Reset();
};
