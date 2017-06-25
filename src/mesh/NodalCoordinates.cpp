
#include "NodalCoordinates.hpp"

#include "PreprocessorExceptions.hpp"

#include <json/json.h>

namespace neon
{
NodalCoordinates::NodalCoordinates(Vector coordinates) : X(coordinates) {}

NodalCoordinates::NodalCoordinates(Json::Value const& mesh_file)
{
    if (mesh_file["Nodes"].empty()) throw EmptyFieldException("Nodes");

    auto const& coordinates = mesh_file["Nodes"][0]["Coordinates"];

    auto const nodes = coordinates.size();

    X.resize(3 * nodes);

    for (auto node = 0; node < nodes; ++node)
    {
        for (auto i = 0; i < 3; ++i)
        {
            X(node * 3 + i) = coordinates[node][i].asDouble();
        }
    }
}

Vector NodalCoordinates::operator[](List const& local_node_list) const
{
    Vector local_coordinates(local_node_list.size() * 3);

    for (auto i = 0; i < local_node_list.size(); ++i)
    {
        local_coordinates.segment<3>(i * 3) = X.segment<3>(local_node_list[i] * 3);
    }
    return local_coordinates;
}
}
