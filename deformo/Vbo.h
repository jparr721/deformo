#pragma once

#include "OpenGL.h"

#include <Eigen/Dense>
#include <string>

template <typename Derived>
void BindVertexAttributeArray(const unsigned int program_id, const std::string& name,
                         const unsigned int buffer, const unsigned int stride,
                         const Eigen::PlainObjectBase<Derived>& data,
                         const bool refresh = true) {
    const unsigned int handle = glGetAttribLocation(program_id, name.c_str());
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    if (refresh) {
        glBufferData(GL_ARRAY_BUFFER, sizeof(Derived::Scalar) * data.size(), data.data(),
                     GL_DYNAMIC_DRAW);
    }
    glVertexAttribPointer(handle, stride, GL_FLOAT, GL_FALSE, 3 * sizeof(Derived::Scalar),
                          nullptr);
    glEnableVertexAttribArray(handle);
}
