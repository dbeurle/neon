
#include "nodal_coordinates.hpp"

#include "exceptions.hpp"
#include "io/json.hpp"

namespace neon
{
nodal_coordinates::nodal_coordinates(matrix3x coordinates) : X(coordinates) {}

nodal_coordinates::nodal_coordinates(json const& mesh_file)
{
    if (mesh_file["Nodes"].empty())
    {
        throw std::domain_error("The mesh file is missing the \"Nodes\" field");
    }

    auto const& coordinates = mesh_file["Nodes"][0]["Coordinates"];

    auto const nodes = coordinates.size();

    X.resize(3, nodes);

    for (std::size_t node{0}; node < nodes; ++node)
    {
        for (auto i = 0; i < 3; ++i)
        {
            X(i, node) = coordinates[node][i];
        }
    }
}
}
