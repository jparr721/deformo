#pragma once

#include <GL/glew.h>
#include <GL/GL.h>
#include <iostream>

inline auto LogErrors(const char* fn) -> void {
    GLenum err;
    for (;;) {
        err = glGetError();

        if (err == GL_NO_ERROR) {
            break;
        }

        std::cerr << "Error in fn: " << fn << ": " << gluErrorString(err) << std::endl;
    }
}
