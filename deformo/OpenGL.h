#pragma once

#ifdef _WIN32
#include <windows.h>
#include <gl/gl.h>
#include <gl/glu.h>
#elif __APPLE__
#define GLEW_STATIC
#include <GL/glew.h>
#endif
