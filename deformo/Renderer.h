#pragma once

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
    unsigned int vao = 0;
    unsigned int vbo = 0;
    unsigned int c_vbo = 0;
    unsigned int ibo = 0;

    GLenum render_mode = GL_FILL;

    Renderer(std::shared_ptr<Mesh> mesh, std::shared_ptr<ShaderProgram> shader_program);
    ~Renderer();

    auto Render() -> void;
    auto SetPositionDisplacement(const Eigen::VectorXf& positions) -> void;
    auto SetColors(const Eigen::VectorXf& colors) -> void;

  private:
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<ShaderProgram> shader_program_;

    std::bitset<4> dirty_;

    auto BuildBuffers() -> void;
    auto ReloadRenderMode() -> void;
    auto ReloadVertexBuffers() -> void;
};
