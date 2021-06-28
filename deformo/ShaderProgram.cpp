#include "ShaderProgram.h"
#include "Utils.h"

ShaderProgram::ShaderProgram() { id = glCreateProgram(); }

auto ShaderProgram::AddShader(GLenum shader_type, const std::string& path)
    -> void {
    const auto shader = std::make_shared<Shader>(shader_type, path);
    assert(shader->Compile() && "FAILED TO COMPILE SHADER");

    shaders.insert({shader_type, shader});
    glAttachShader(id, shader->id);
    LogErrors("AddShader");
}

auto ShaderProgram::Link() const -> void { glLinkProgram(id); }

ShaderProgram::~ShaderProgram() {
    for (const auto& [k, shader] : shaders) {
        glDeleteShader(shader->id);
    }
    glDeleteProgram(id);
}

auto ShaderProgram::Release() const -> void { glUseProgram(0); }

auto ShaderProgram::SetMatrixUniformIdentity() -> void { glLoadIdentity(); }

auto ShaderProgram::Bind() const -> void { glUseProgram(id); }

auto ShaderProgram::SetMatrixUniform(int location,
                                     const Eigen::MatrixXf& uniform) -> void {
    glUniformMatrix4fv(location, 1, GL_FALSE, uniform.data());
}

auto ShaderProgram::UniformLocation(const std::string& name) -> int {
    return glGetUniformLocation(id, name.c_str());
}
