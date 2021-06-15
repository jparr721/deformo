#pragma once

#include "OpenGL.h"
#include <string>
#include <utility>

class Shader {
  public:
    GLuint id = -1;
    GLenum type;

    Shader(const GLenum type, std::string shader_path)
        : type(type), shader_path_(std::move(shader_path)) {}

    auto Compile() -> bool;

  private:
    const std::string shader_path_;
};
