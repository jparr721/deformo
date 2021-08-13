#pragma once

#include "OpenGL.h"

#include <Eigen/Dense>

template <typename Derived>
void BindElementArrayObject(unsigned int& buffer,
                            const Eigen::PlainObjectBase<Derived>& data,
                            bool reload = false) {
    if (!reload) {
        glGenBuffers(1, &buffer);
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Derived::Scalar) * data.size(),
                 data.data(), GL_STATIC_DRAW);
}
