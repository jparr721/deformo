#include "MarchingCubes.h"
#include <map>

namespace {
struct Position {
    unsigned int index;
    Vector3r position;
};

auto Interpolate(const Real x1, const Real y1, const Real z1, const Real x2,
                 const Real y2, const Real z2, const Real val1, const Real val2,
                 const Real iso_level) -> Vector3r {
    const auto mu_interpolate = [](const Real v1, const Real v2,
                                   const Real mu) -> Real {
        return v1 + mu * (v2 - v1);
    };

    const Real mu = (iso_level - val1) / (val1 - val2);
    return Vector3r(mu_interpolate(x1, x2, mu), mu_interpolate(y1, y2, mu),
                    mu_interpolate(z1, z2, mu));
}

auto Intersection(const unsigned int x, const unsigned int y,
                  const unsigned int z, const unsigned int edge_number,
                  const unsigned int cell_length, const Real iso_level,
                  const Tensor3r& scalar_field) -> Vector3r {
    unsigned int v1x = x;
    unsigned int v1y = y;
    unsigned int v1z = z;
    unsigned int v2x = x;
    unsigned int v2y = y;
    unsigned int v2z = z;

    switch (edge_number) {
    case 0:
        v2y += 1;
        break;
    case 1:
        v1y += 1;
        v2x += 1;
        v2y += 1;
        break;
    case 2:
        v1x += 1;
        v1y += 1;
        v2x += 1;
        break;
    case 3:
        v1x += 1;
        break;
    case 4:
        v1z += 1;
        v2y += 1;
        v2z += 1;
        break;
    case 5:
        v1y += 1;
        v1z += 1;
        v2x += 1;
        v2y += 1;
        v2z += 1;
        break;
    case 6:
        v1x += 1;
        v1y += 1;
        v1z += 1;
        v2x += 1;
        v2z += 1;
        break;
    case 7:
        v1x += 1;
        v1z += 1;
        v2z += 1;
        break;
    case 8:
        v2z += 1;
        break;
    case 9:
        v1y += 1;
        v2y += 1;
        v2z += 1;
        break;
    case 10:
        v1x += 1;
        v1y += 1;
        v2x += 1;
        v2y += 1;
        v2z += 1;
        break;
    case 11:
        v1x += 1;
        v2x += 1;
        v2z += 1;
        break;
    default:
        assert(edge_number <= 11 && edge_number >= 0 &&
               "INVALID EDGE NUMBER FOUND");
    }

    const Real x1 = v1x * cell_length;
    const Real y1 = v1y * cell_length;
    const Real z1 = v1z * cell_length;
    const Real x2 = v2x * cell_length;
    const Real y2 = v2y * cell_length;
    const Real z2 = v2z * cell_length;

    const Real val1 = scalar_field.At(v1z, v1y, v1x);
    const Real val2 = scalar_field.At(v2z, v2y, v2x);
    return Interpolate(x1, y1, z1, x2, y2, z2, val1, val2, iso_level);
}

auto VertexId(const unsigned int x, const unsigned int y, const unsigned int z,
              const unsigned int n_cells_x, const unsigned int n_cells_y)
    -> unsigned int {
    return 3 *
           (z * (n_cells_y + 1) * (n_cells_x + 1) + y * (n_cells_x + 1) + x);
}

auto EdgeId(const unsigned int x, const unsigned int y, const unsigned int z,
            const unsigned int n_cells_x, const unsigned int n_cells_y,
            const unsigned int edge_number) -> unsigned int {
    switch (edge_number) {
    case 0:
        return VertexId(x, y, z, n_cells_x, n_cells_y) + 1;
    case 1:
        return VertexId(x, y + 1, z, n_cells_x, n_cells_y);
    case 2:
        return VertexId(x + 1, y, z, n_cells_x, n_cells_y) + 1;
    case 3:
        return VertexId(x, y, z, n_cells_x, n_cells_y);
    case 4:
        return VertexId(x, y, z + 1, n_cells_x, n_cells_y) + 1;
    case 5:
        return VertexId(x, y + 1, z + 1, n_cells_x, n_cells_y);
    case 6:
        return VertexId(x + 1, y, z + 1, n_cells_x, n_cells_y) + 1;
    case 7:
        return VertexId(x, y, z + 1, n_cells_x, n_cells_y);
    case 8:
        return VertexId(x, y, z, n_cells_x, n_cells_y) + 2;
    case 9:
        return VertexId(x, y + 1, z, n_cells_x, n_cells_y) + 2;
    case 10:
        return VertexId(x + 1, y + 1, z, n_cells_x, n_cells_y) + 2;
    case 11:
        return VertexId(x + 1, y, z, n_cells_x, n_cells_y) + 2;
    default:
        return -1;
    }
}

void CalculateFacesAndVertices(MatrixXr& V, Eigen::MatrixXi& F,
                               std::vector<Eigen::Vector3i>& triangles,
                               std::map<unsigned int, Position>& vertices) {
    unsigned int index = 0;
    for (auto& [_, position] : vertices) {
        position.index = index;
        ++index;
    }

    for (auto& triangle : triangles) {
        const unsigned int new_index = vertices.at(triangle.x()).index;
        triangle.x() = new_index;
        triangle.y() = new_index;
        triangle.z() = new_index;
    }

    V.resize(vertices.size(), 3);
    unsigned int i = 0;
    for (const auto& [_, position] : vertices) {
        V.row(i) = position.position;
        ++i;
    }

    F.resize(triangles.size(), 3);
    i = 0;
    for (const auto& triangle : triangles) {
        F.row(i) = triangle;
        ++i;
    }
}
} // namespace

void GenerateFacesFromScalarField(MatrixXr& V, Eigen::MatrixXi& F,
                                  const Real iso_level, const Real cell_length,
                                  const Tensor3r& scalar_field) {
    std::vector<Eigen::Vector3i> triangles;
    std::map<unsigned int, Position> vertices;

    const unsigned int n_cells_z = scalar_field.Dimension(0);
    const unsigned int n_cells_y = scalar_field.Dimension(1);
    const unsigned int n_cells_x = scalar_field.Dimension(2);

    const auto insert_at_edge = [&](const unsigned int x, const unsigned int y,
                                    const unsigned int z,
                                    const unsigned int edge_number) -> void {
        const Vector3r point = Intersection(x, y, z, edge_number, cell_length,
                                            iso_level, scalar_field);
        const unsigned int id =
            EdgeId(x, y, z, n_cells_x, n_cells_y, edge_number);
        vertices.insert({id, Position{0, point}});
    };

    for (unsigned int z = 0; z < n_cells_z; ++z) {
        for (unsigned int y = 0; y < n_cells_y; ++y) {
            for (unsigned int x = 0; x < n_cells_x; ++x) {
                // Calculate lookup index in table
                unsigned int table_index = 0;
                if (scalar_field(x, y, z) < iso_level) {
                    table_index |= 1;
                }

                if (scalar_field(x, y + 1, z) < iso_level) {
                    table_index |= 2;
                }

                if (scalar_field(x + 1, y + 1, z) < iso_level) {
                    table_index |= 4;
                }

                if (scalar_field(x + 1, y, z) < iso_level) {
                    table_index |= 8;
                }

                if (scalar_field(x, y, z + 1) < iso_level) {
                    table_index |= 16;
                }

                if (scalar_field(x, y + 1, z + 1) < iso_level) {
                    table_index |= 32;
                }

                if (scalar_field(x + 1, y + 1, z + 1) < iso_level) {
                    table_index |= 64;
                }

                if (scalar_field(x + 1, y, z + 1) < iso_level) {
                    table_index |= 128;
                }

                // Triangulate in this cell
                if (edge_table.at(table_index) != 0) {
                    if (edge_table.at(table_index) & 8) {
                        constexpr unsigned int edge_number = 3;
                        insert_at_edge(x, y, z, edge_number);
                    }

                    if (edge_table.at(table_index) & 1) {
                        constexpr unsigned int edge_number = 0;
                        insert_at_edge(x, y, z, edge_number);
                    }

                    if (edge_table.at(table_index) & 256) {
                        constexpr unsigned int edge_number = 8;
                        insert_at_edge(x, y, z, edge_number);
                    }

                    if (x == n_cells_x - 1) {
                        if (edge_table.at(table_index) & 4) {
                            constexpr unsigned int edge_number = 2;
                            insert_at_edge(x, y, z, edge_number);
                        }

                        if (edge_table.at(table_index) & 2048) {
                            constexpr unsigned int edge_number = 11;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    if (y == n_cells_y - 1) {
                        if (edge_table.at(table_index) & 2) {
                            constexpr unsigned int edge_number = 1;
                            insert_at_edge(x, y, z, edge_number);
                        }

                        if (edge_table.at(table_index) & 512) {
                            constexpr unsigned int edge_number = 9;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    if (z == n_cells_z - 1) {
                        if (edge_table.at(table_index) & 16) {
                            constexpr unsigned int edge_number = 4;
                            insert_at_edge(x, y, z, edge_number);
                        }

                        if (edge_table.at(table_index) & 128) {
                            constexpr unsigned int edge_number = 7;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    if ((x == n_cells_x - 1) && (y == n_cells_y - 1)) {
                        if (edge_table.at(table_index) & 1024) {
                            constexpr unsigned int edge_number = 10;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    if ((x == n_cells_x - 1) && (z == n_cells_z - 1)) {
                        if (edge_table.at(table_index) & 64) {
                            constexpr unsigned int edge_number = 6;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    if ((y == n_cells_y - 1) && (z == n_cells_z - 1)) {
                        if (edge_table.at(table_index) & 32) {
                            constexpr unsigned int edge_number = 5;
                            insert_at_edge(x, y, z, edge_number);
                        }
                    }

                    for (unsigned int i = 0;
                         triangle_table.at(table_index).at(i) != -1; i += 3) {
                        Eigen::Vector3i triangle;
                        triangle.x() =
                            EdgeId(x, y, z, n_cells_x, n_cells_y,
                                   triangle_table.at(table_index).at(i));
                        triangle.y() =
                            EdgeId(x, y, z, n_cells_x, n_cells_y,
                                   triangle_table.at(table_index).at(i + 1));
                        triangle.z() =
                            EdgeId(x, y, z, n_cells_x, n_cells_y,
                                   triangle_table.at(table_index).at(i + 2));
                        triangles.emplace_back(triangle);
                    }
                }
            }
        }
    }

    assert(triangles.empty() && "NO RECOREDED TRIANGLES");
    assert(vertices.empty() && "NO RECORDED VERTICES");
    CalculateFacesAndVertices(V, F, triangles, vertices);
}
