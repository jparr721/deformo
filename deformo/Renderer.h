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
    int mvp = -1;

    Renderer(std::shared_ptr<Mesh> mesh,
             std::shared_ptr<ShaderProgram> shader_program,
             std::shared_ptr<Camera<Real>> camera);
    ~Renderer();

    auto Render() -> void;
    auto Resize(int width, int height) const -> void;
    auto SetPositionDisplacement(const VectorXr& positions) -> void;
    auto SetColors(const VectorXr& colors) -> void;
    auto SetRenderMode(GLenum mode) -> void;
    void SetTetgenFlags(const std::string& flags);
    void SetCutPlane(Real cut_plane);
    void SetCutPlaneAxis(SliceAxis cut_plane_axis);

  private:
    GLenum render_mode_ = GL_LINE;

    std::bitset<4> dirty_;

    std::shared_ptr<Camera<Real>> camera_;
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<ShaderProgram> shader_program_;

    auto BuildBuffers() -> void;
    auto ReloadRenderMode() -> void;
    auto ReloadVertexBuffers() -> void;
};
