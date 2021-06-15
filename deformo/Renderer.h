#pragma once

#include "Mesh.h"
#include "OpenGL.h"
//#include "Shader.h"
//#include <bitset>
//#include <memory>
//
//enum DirtyStatus {
//    POSITIONS = 0,
//    COLORS,
//    TEXTURES,
//    RENDER_MODE,
//};
//
//class Renderer {
//  public:
//    unsigned int vao = 0;
//    unsigned int vbo = 0;
//    unsigned int c_vbo = 0;
//    unsigned int ibo = 0;
//
//    GLenum render_mode = GL_FILL;
//
//    Renderer(std::shared_ptr<Mesh> mesh, std::shared_ptr<Shader> shader);
//    ~Renderer();
//
//    auto Render() -> void;
//    auto SetDirty(DirtyStatus status) -> void;
//
//  private:
//    std::shared_ptr<Mesh> mesh_;
//    std::shared_ptr<Shader> shader_;
//
//    std::bitset<4> dirty_;
//
//    auto ReloadRenderMode() -> void;
//    auto ReloadVertexBuffers() -> void;
//};
