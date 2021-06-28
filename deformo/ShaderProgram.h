#pragma once

#include <QMatrix4x4>
#include "OpenGL.h"
#include "Shader.h"
#include <Eigen/Dense>
#include <string>
#include <unordered_map>

class ShaderProgram {

  public:
    int id = -1;

    std::unordered_map<GLenum, std::shared_ptr<Shader>> shaders;

    ShaderProgram();
    ~ShaderProgram();

    auto AddShader(GLenum shader_type, const std::string& path) -> void;
    auto Link() const -> void;
    auto Bind() const -> void;
    auto Release() const -> void;

    auto SetMatrixUniformIdentity() -> void;
    auto SetMatrixUniform(int location, const Eigen::MatrixXf& uniform) -> void;

    auto UniformLocation(const std::string& name) -> int;
};
