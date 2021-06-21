#pragma once

#include "Camera.h"
#include "Mesh.h"
#include "OpenGL.h"
#include "ShaderProgram.h"
#include <bitset>
#include <memory>

enum DirtyStatus {
    positions = 0,
    colors,
    textures,
    render_mode,
};

class Renderer {
  public:
    // Geometry Rendering System
    unsigned int vao = 0;
    unsigned int vbo = 0;
    unsigned int c_vbo = 0;
    unsigned int ibo = 0;

    // View Coordinate System
    int m = -1;
    int v = -1;
    int p = -1;

    Renderer(std::shared_ptr<Mesh> mesh,
             std::shared_ptr<ShaderProgram> shader_program,
             std::shared_ptr<Camera> camera);
    ~Renderer();

    auto Render() -> void;
    auto SetPositionDisplacement(const Eigen::VectorXf& positions) -> void;
    auto SetColors(const Eigen::VectorXf& colors) -> void;
    auto SetRenderMode() -> void;
    void SetTetgenFlags(const std::string& flags);
    void SetCutPlane(float cut_plane);
    void SetCutPlaneAxis(CutPlaneAxis cut_plane_axis);

  private:
    GLenum render_mode_ = GL_LINE;

    std::bitset<4> dirty_;

    std::shared_ptr<Camera> camera_;
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<ShaderProgram> shader_program_;

    auto BuildBuffers() -> void;
    auto ReloadRenderMode() -> void;
    auto ReloadVertexBuffers() -> void;
};
