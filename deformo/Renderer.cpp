#include "Renderer.h"

#include "Configuration.h"
#include "Vao.h"
#include "Vbo.h"

#include <utility>

Renderer::Renderer(std::shared_ptr<Mesh> mesh,
                   std::shared_ptr<ShaderProgram> shader_program)
    : mesh_(std::move(mesh)), shader_program_(std::move(shader_program)) {
    shader_program_->AddShader(
        GL_VERTEX_SHADER,
        env_config->Value(config_keys::kVertexShaderLocation));
    shader_program_->AddShader(
        GL_FRAGMENT_SHADER,
        env_config->Value(config_keys::kFragmentShaderLocation));

    shader_program_->Link();
    shader_program_->Bind();
    BuildBuffers();
    shader_program_->Release();
}

Renderer::~Renderer() {
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &c_vbo);
    glDeleteBuffers(1, &ibo);
    glDeleteVertexArrays(1, &vao);
}

auto Renderer::Render() -> void {
    ReloadRenderMode();
    ReloadVertexBuffers();
    shader_program_->Bind();

    glDrawElements(GL_TRIANGLES, mesh_->faces.size(), GL_UNSIGNED_INT, nullptr);

    shader_program_->Release();
}

auto Renderer::SetPositionDisplacement(const Eigen::VectorXf& positions)
    -> void {
    mesh_->Update(positions);
    dirty_.set(DirtyStatus::positions);
}

auto Renderer::SetColors(const Eigen::VectorXf& colors) -> void {
    mesh_->colors = colors;
    dirty_.set(DirtyStatus::colors);
}

auto Renderer::BuildBuffers() -> void {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &c_vbo);

    BindVertexAttributeArray(shader_program_->id, "position", vbo, 3,
                             mesh_->positions);
    BindVertexAttributeArray(shader_program_->id, "colors", c_vbo, 3,
                             mesh_->colors);
    BindElementArrayObject(ibo, mesh_->faces);
}

auto Renderer::ReloadRenderMode() -> void {
    if (dirty_[DirtyStatus::render_mode]) {
        glPolygonMode(GL_FRONT_AND_BACK, render_mode);
    }
}

auto Renderer::ReloadVertexBuffers() -> void {
    if (dirty_[DirtyStatus::positions]) {
        BindVertexAttributeArray(shader_program_->id, "position", vbo, 3,
                                 mesh_->positions);
        dirty_.flip(DirtyStatus::positions);
    }
    if (dirty_[DirtyStatus::colors]) {
        BindVertexAttributeArray(shader_program_->id, "colors", c_vbo, 3,
                                 mesh_->colors);
        dirty_.flip(DirtyStatus::colors);
    }
}
