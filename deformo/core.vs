#version 330

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 color;

out vec4 v_color;

uniform mat4 m, v, p;

void main() {
    v_color = vec4(color, 1.0);
    gl_Position = p * vec4(position, 1.0);
}