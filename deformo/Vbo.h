#pragma once

#include "OpenGL.h"
#include <Eigen/Dense>

template <typename Derived> class Vbo {
  public:
    using Data = Eigen::MatrixBase<Derived>;
    bool dirty = false;
    const GLuint handle;
    const GLenum vbo_type;
    const GLenum usage;

    Data data;

    Vbo(GLuint handle, GLenum vbo_type, GLenum usage)
        : handle(handle), vbo_type(vbo_type), usage(usage) {}

    auto Update(const Data& input) -> void {
        data += input;
        dirty = true;
    }

    auto Refresh() -> void {
        if (!dirty) {
            return;
        }

        BindVertexAttributeArray(data);

        dirty = false;
    }

    auto BindVertexAttributeArray(const Data& data) -> void {
        glBindBuffer(vbo_type, handle);
    }
};
