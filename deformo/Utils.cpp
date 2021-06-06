#include "Utils.h"

namespace utils {
void FindMaxVertices(std::vector<unsigned>& indices,
                     const Eigen::VectorXf& positions) {
    float max_y = -1e-10f;
    // Iterate over each y value in the position vector
    for (int i = 1; i < positions.rows(); i += 3) {
        max_y = std::fmax(positions(i), max_y);
    }

    for (unsigned int i = 1; i < positions.rows(); i += 3) {
        if (positions(i) >= max_y) {
            // Add the x-oriented value so we can fetch the whole group later
            indices.push_back((i - 1));
        }
    }
}
} // namespace utils