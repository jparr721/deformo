#include "Renderer.h"

#include "Dotenv.h"
#include "Vao.h"
#include "Vbo.h"

#include <utility>

Renderer::Renderer(std::shared_ptr<Mesh> mesh,
                   std::shared_ptr<ShaderProgram> shader_program,
                   std::shared_ptr<Camera<Real>> camera)
    : camera_(std::move(camera)), mesh_(std::move(mesh)),
      shader_program_(std::move(shader_program)) {
    shader_program_->AddShader(
        GL_VERTEX_SHADER,
        env_config->Value(config_keys::kVertexShaderLocation));
    shader_program_->AddShader(
        GL_FRAGMENT_SHADER,
        env_config->Value(config_keys::kFragmentShaderLocation));

    shader_program_->Link();
    shader_program_->Bind();

    mvp = shader_program_->UniformLocation("mvp");

    BuildBuffers();
    shader_program_->Release();

    // Set render mode so it runs on first render
    dirty_.set(DirtyStatus::render_mode);
}

Renderer::~Renderer() {
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &c_vbo);
    glDeleteBuffers(1, &ibo);
    glDeleteVertexArrays(1, &vao);
}

auto Renderer::Render() -> void {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ReloadRenderMode();
    ReloadVertexBuffers();

    shader_program_->Bind();

    shader_program_->SetMatrixUniformIdentity();
    shader_program_->SetMatrixUniform(mvp, camera_->ToViewMatrix());

    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, mesh_->FacesSize(), GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);

    shader_program_->Release();

    LogErrors("Renderer::Render");
}

auto Renderer::Resize(int width, int height) -> void {
    glViewport(0, 0, width, height);
    shader_program_->Bind();
    camera_->Resize(width, height);
    shader_program_->SetMatrixUniform(mvp, camera_->GetProjectionMatrix());
    shader_program_->Release();
    LogErrors("Renderer::Resize");
}

auto Renderer::SetPositionDisplacement(const VectorXr& positions) -> void {
    mesh_->Update(positions);
    dirty_.set(DirtyStatus::positions);
}

auto Renderer::SetColors(const VectorXr& colors) -> void {
    mesh_->colors = colors;
    dirty_.set(DirtyStatus::colors);
}

auto Renderer::SetRenderMode(const GLenum mode) -> void {
    render_mode_ = mode;
    dirty_.set(DirtyStatus::render_mode);
}

void Renderer::SetTetgenFlags(const std::string& flags) {}

void Renderer::SetCutPlane(Real cut_plane) { mesh_->SetSliceValue(cut_plane); }

void Renderer::SetCutPlaneAxis(SliceAxis cut_plane_axis) {
    mesh_->SetSliceAxis(cut_plane_axis);
}

auto Renderer::BuildBuffers() -> void {
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &c_vbo);

    BindVertexAttributeArray(shader_program_->id, "position", vbo, 3,
                             mesh_->positions);
    BindVertexAttributeArray(shader_program_->id, "color", c_vbo, 4,
                             mesh_->colors);
    BindElementArrayObject(ibo, mesh_->faces);
    LogErrors("Renderer::BuildBuffers");
}

auto Renderer::ReloadRenderMode() -> void {
    if (dirty_[DirtyStatus::render_mode]) {
        glPolygonMode(GL_FRONT_AND_BACK, render_mode_);
        dirty_.flip(DirtyStatus::render_mode);
    }
}

auto Renderer::ReloadVertexBuffers() -> void {
    BindVertexAttributeArray(shader_program_->id, "position", vbo, 3,
                             mesh_->positions);
    // if (dirty_[DirtyStatus::colors]) {
    BindVertexAttributeArray(shader_program_->id, "color", c_vbo, 4,
                             mesh_->colors);
    //    dirty_.flip(DirtyStatus::colors);
    //}
}
