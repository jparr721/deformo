attribute highp vec4 position;
attribute lowp vec4 color;

varying lowp vec4 v_color;

uniform highp mat4 projection;
uniform highp mat4 model;
uniform highp mat4 view;

void main() {
    v_color = color;
    gl_Position = projection * position;
}