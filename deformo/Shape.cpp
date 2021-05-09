#include "Shape.h"

void MarchCubes(double surface_value, const std::vector<double>& scalar_field,
                std::vector<Eigen::Vector3d>& triangles,
                std::vector<Eigen::Vector3d>& normals) {
  const double step_size = 1.0 / scalar_field.size();
  for (int x = 0; x < scalar_field.size(); ++x) {
    for (int y = 0; y < scalar_field.size(); ++y) {
      for (int z = 0; z < scalar_field.size(); ++z) {
        MarchCube, (x * step_size, y * step_size, z * step_size, step_size, , );
      }
    }
  }
}

void MarchCube(double x, double y, double z, double surface_value, double scale,
               std::vector<Eigen::Vector3d>& triangles,
               std::vector<Eigen::Vector3d>& normals) {
  int flag_index;

  float offset;

  std::vector<float> current_vertex_value;

  for (int i = 0; i < shape::kNCubeVertices; ++i) {
  }
}
