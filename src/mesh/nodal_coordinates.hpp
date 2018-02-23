
#pragma once

#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"

#include "io/json.hpp"

namespace neon
{
/** nodal_coordinates.hpp stores the list of coordinates of a discretized geometry */
class nodal_coordinates
{
public:
    nodal_coordinates() = default;

    /** Construct with a list of coordinates */
    explicit nodal_coordinates(matrix3x const coordinates);

    /** Construct with a list of coordinates in json format */
    explicit nodal_coordinates(json const& mesh_file);

    [[nodiscard]] auto size() const { return X.cols(); }

    [[nodiscard]] matrix3x const& coordinates() const { return X; }

    /** @return the coordinates using fancy indexing */
    template <typename indices_type>
    [[nodiscard]] auto coordinates(indices_type const local_node_view) const
    {
        return X(Eigen::placeholders::all, local_node_view);
    }

protected:
    matrix3x X; //!< Reference configuration encoded as (x1, y1, z1, x2, y2, z2)
};

/**
 *
 */
template <typename traits>
class mesh_coordinates
{
    /** fixed size coordinates depending on the mathematical model */
    static auto constexpr fixed_size = traits::size;

    using coordinate_t = Eigen::Matrix<double, fixed_size, Eigen::Dynamic>;

public:
    /** Construct with a list of coordinates */
    explicit mesh_coordinates(coordinate_t const coordinates);

    /** Construct with a list of coordinates in json format */
    explicit mesh_coordinates(json const& mesh_file);

    [[nodiscard]] auto size() const { return X.cols(); }

    [[nodiscard]] coordinate_t const& coordinates() const { return X; }

        /** @return the coordinates using fancy indexing */
        [[nodiscard]] coordinate_t coordinates(index_view local_node_view) const
    {
        return X(Eigen::placeholders::all, local_node_view);
    }

protected:
    coordinate_t X; //!< Reference configuration encoded as (x1, y1, z1, x2, y2, z2)
};

template <typename traits>
mesh_coordinates<traits>::mesh_coordinates(mesh_coordinates::coordinate_t coordinates)
    : X(coordinates)
{
}

template <typename traits>
mesh_coordinates<traits>::mesh_coordinates(json const& mesh_file)
{
    if (mesh_file["Nodes"].is_null())
    {
        throw std::domain_error("The mesh file is missing the \"Nodes\" field");
    }

    auto const& input_coordinates = mesh_file["Nodes"][0]["Coordinates"];

    auto const nodes = input_coordinates.size();

    X.resize(fixed_size, nodes);

    for (std::int64_t node{0}; node < nodes; ++node)
    {
        for (auto i{0}; i < fixed_size; ++i)
        {
            X(i, node) = input_coordinates[node][i];
        }
    }
}
}
