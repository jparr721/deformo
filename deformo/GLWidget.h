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

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

  public:
    unsigned int vao = 0;
    unsigned int vbo = 0;
    unsigned int c_vbo = 0;
    unsigned int ibo = 0;

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

  public slots:
    void Cleanup();
    void Update();

    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void SetSliceAxis(const QString& value);
    void SetSliceValue(float value);
    void SetNodalMass(float value);
    void SetPoissonsRatio(double value);
    void SetYoungsModulus(double value);
    void SetTimestepSize(double value);

    // Simulation Settings -- Damping Parameters
    void SetRayleighLambda(double value);
    void SetRayleighMu(double value);

    // Simulation Settings -- Run Simulation Button
    void RunSimulationButtonPressed();

    // Render Properties Window -- Tetgen Flags
    void SetTetgenFlags(const QString& value);

    // Render Properties Window -- Re Render Button
    void RenderSimulationButtonPressed();

  signals:
    // Simulation Settings Window
    // Simulation Settings -- FEA Parameters
    void OnSliceAxisChange(const QString& value);
    void OnSliceValueChange(float value);
    void OnNodalMassChange(float value);
    void OnPoissonsRatioChange(double value);
    void OnYoungsModulusChange(double value);
    void OnTimestepSizeChange(double value);

    // Simulation Settings -- Damping Parameters
    void OnRayleighLambdaChange(double value);
    void OnRayleighMuChange(double value);

    // Render Properties Window
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
    std::shared_ptr<Camera> camera_;
    bool simulating_ = false;

    void BuildMesh(float cut_plane = Mesh::kNoCutPlane);
    void BuildPhysicsEngine();

    void Reset();
};
