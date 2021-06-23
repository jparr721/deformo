#include "Utils.h"
#include <sstream>

namespace utils {
std::string QStringToString(const QString& string) {
    return string.toUtf8().constData();
}

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
auto OpenFile(const std::string& filename) -> std::ifstream {
    std::ifstream input;
    input.open(filename);
    assert(input.good() && "FAILED TO READ FILE");

    return input;
}

auto ReadFile(const std::string& path) -> std::string {
    std::ostringstream sstr;
    const auto stream = OpenFile(path);
    sstr << stream.rdbuf();
    return sstr.str();
}

auto Split(const std::string& input, const std::string& delimiter)
    -> std::pair<std::string, std::string> {
    const auto delimiter_index = input.find(delimiter);
    assert(delimiter_index != input.npos && "MALFORMED CONFIG");
    const std::string key = input.substr(0, delimiter_index);
    const std::string value = input.substr(delimiter_index + 1);
    return std::make_pair(key, value);
}
} // namespace utils