#pragma once

#include "ShaderProgram.h"
#include "Vao.h"
#include "Vbo.h"

#include <Eigen/Dense>
#include <QKeyEvent>
#include <QOpenGLFunctions>
#include <QOpenGLWidget>
#include <QTimer>
#include <memory>

#include "Camera.h"
#include "Input.h"
#include "Mesh.h"
#include "Renderer.h"
#include "Simulation.h"

class WindowController;

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

  public:
    std::shared_ptr<ShaderProgram> shader_program;

    // Toggleable Wire Mesh
    GLenum render_style = GL_LINE;

    // Mesh object for 3D geometry
    std::shared_ptr<Mesh> mesh;

    // Renderer for loading the sim geometry
    std::unique_ptr<Renderer> renderer;

    // LinearTetrahedral Object
    std::unique_ptr<Simulation> sim;

    explicit GLWidget(QWidget* parent = nullptr);
    ~GLWidget() override;

    void resizeGL(int width, int height) override;
    void SetController(const std::shared_ptr<WindowController> controller) {
        controller_ = controller;
    }

  public slots:
    void Update();

  protected:
    void initializeGL() override;
    void paintGL() override;

    // Mouse Shenanigans
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

  private:
    QTimer* draw_timer_;
    std::unique_ptr<Input> input_;
    std::shared_ptr<Camera<Real>> camera_;
    std::shared_ptr<WindowController> controller_;
    bool simulating_ = false;
};
