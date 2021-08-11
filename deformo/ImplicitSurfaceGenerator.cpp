#include "ImplicitSurfaceGenerator.h"

_Cube::_Cube(int rows, int cols, int layers,
             const Vector3<unsigned int>& starting_index)
    : rows_(rows), cols_(cols), layers_(layers) {
    const unsigned int c_x = starting_index.x();
    const unsigned int c_y = starting_index.y();
    const unsigned int c_z = starting_index.z();

    for (unsigned int layer = 0; layer < layers_; ++layer) {
        for (unsigned int col = 0; col < cols_; ++col) {
            for (unsigned int row = 0; row < rows_; ++row) {
                indices_.emplace_back(c_x + row, c_y + col, c_z + layer);
            }
        }
    }

    // Remove any duplicates
    std::sort(indices_.begin(), indices_.end(), UnsignedVectorComparison);
    indices_.erase(
        std::unique(indices_.begin(), indices_.end(), UnsignedVectorEqual),
        indices_.end());
}

bool _Cube::operator==(const _Cube& rhs) const {
    for (const Vector3<unsigned int>& lhs_index : indices_) {
        for (const Vector3<unsigned int>& rhs_index : rhs.Indices()) {
            if (lhs_index == rhs_index) {
                return true;
            }
        }
    }

    return false;
}
