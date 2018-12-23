
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include "io/json.hpp"

namespace neon
{
namespace v2
{
/// \tparam N dimension of coordinate space
template <int N>
class nodal_coordinates
{
public:
    /// Fixed size coordinates depending on the mathematical model
    static auto constexpr size_type = N;
    /// Coordinate storage type
    using value_type = matrixdx<N>;

public:
    /// Construct with a list of coordinates
    explicit nodal_coordinates(value_type const coordinates);

    /// Construct with a list of coordinates in json format
    explicit nodal_coordinates(json const& mesh_file);

    /// \return number of nodal coordinates
    [[nodiscard]] auto size() const noexcept { return X.cols(); }

    /// \return constant access to all coordinates
    [[nodiscard]] value_type const& coordinates() const noexcept { return X; }

    /// \return the coordinates using an index view or single index
    template <typename indices_type>
    [[nodiscard]] auto coordinates(indices_type const local_node_view) const noexcept
    {
        return X(Eigen::all, local_node_view);
    }

protected:
    /// Reference configuration
    value_type X;
};

template <int size>
nodal_coordinates<size>::nodal_coordinates(nodal_coordinates::value_type coordinates)
    : X(coordinates)
{
}

template <int N>
nodal_coordinates<N>::nodal_coordinates(json const& mesh_file)
{
    if (mesh_file["Nodes"].is_null())
    {
        throw std::domain_error("The mesh file is missing the \"Nodes\" field");
    }

    auto const& input_coordinates = mesh_file["Nodes"][0]["Coordinates"];

    auto const nodes = input_coordinates.size();

    X.resize(N, nodes);

    for (std::size_t node{0}; node < nodes; ++node)
    {
        for (auto i{0}; i < N; ++i)
        {
            X(i, node) = input_coordinates[node][i];
        }
    }
}

extern template class nodal_coordinates<2>;
extern template class nodal_coordinates<3>;
}
}
