attribute highp vec4 position_attribute;
attribute lowp vec4 color_attribute;

varying lowp vec4 color;
uniform highp mat4 projection_matrix;

void main() {
    color = color_attribute;
    gl_Position = projection_matrix * position_attribute;
}