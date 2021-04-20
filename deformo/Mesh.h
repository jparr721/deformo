#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

struct Vertex {
  Eigen::Vector3d position;
  Eigen::Vector3d normal;
  Eigen::Vector2d texture_coordinates;
};

struct Texture {
  // ID of the texture in the file
  unsigned int id;

  // Diffuse, specular, etc...
  std::string type;
};

class Mesh {
 public:
  std::vector<Vertex> vertices;
  std::vector<unsigned int> indices;
  std::vector<Texture> textures;

  Mesh(std::vector<Vertex> vertices_, std::vector<unsigned int> indices_, std::vector<Texture> textures_);

  void Draw();

 private:
  unsigned int vao = 0;
  unsigned int vbo = 0;
  unsigned int ebo = 0;

  void SetupMesh();
};
