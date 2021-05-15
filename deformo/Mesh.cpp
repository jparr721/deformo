#include "Mesh.h"

#include <Eigen/Core>

Mesh::Mesh(const Eigen::MatrixXf& V,
 const Eigen::MatrixXf& F,
           const Eigen::MatrixXf& T)
    : faces(F), tetrahedrals(T) {
  InitializeVertexPositions(V);
  vbo = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(vertices.data(), size_bytes());
  vbo.release();

  ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
  ibo.create();
  ibo.bind();
  // Indices won't be changing
  ibo.setUsagePattern(QOpenGLBuffer::StaticDraw);
  ibo.release();
};

void Mesh::InitializeVertexPositions(const Eigen::MatrixXf& V) {
  vertices.resize(V.rows() * V.cols());
  Eigen::MatrixXf Vt = V.transpose();
  vertices = Eigen::Map<Eigen::VectorXf>(Vt.data(), Vt.rows() * Vt.cols());

  colors.resize(vertices.rows());
  for (int i = 0; i < vertices.rows(); i += 3) {
    colors.segment(i, 3) << 1.0f, 0.f, 0.f;
  }
}

void Mesh::Render() {
  // Add updated vertex coordinates
  vbo.bind();
  vbo.write(0, vertices.data(), size_bytes());
  vbo.release();

  // Render
  vao.bind();
  glDrawArrays(GL_TRIANGLES, 0, size());
  vao.release();
}
