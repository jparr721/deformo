#include "Shader.h"
#include "Utils.h"

auto Shader::Compile() -> bool {
    const std::string src = utils::ReadFile(shader_path_);
    const char* shader_source = src.c_str();
    id = glCreateShader(type);
    glShaderSource(id, 1, &shader_source, nullptr);
    glCompileShader(id);

    GLint is_linked;
    glGetShaderiv(id, GL_COMPILE_STATUS, &is_linked);
    if (is_linked == GL_FALSE) {
        GLsizei log_length = 0;
        GLchar message[1024];
        glGetShaderInfoLog(id, 1024, &log_length, message);
        std::cerr << message << std::endl;
        return false;
    }

    LogErrors("Compile");

    return is_linked == GL_TRUE;
}
