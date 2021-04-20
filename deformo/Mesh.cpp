#include "Mesh.h"

Mesh::Mesh(std::vector<Vertex> vertices_, std::vector<unsigned int> indices_,
           std::vector<Texture> textures_) :vertices(vertices_), indices(indices_), textures(textures_) {
  SetupMesh();
}

void Mesh::Draw() {}

void Mesh::SetupMesh() {
}
