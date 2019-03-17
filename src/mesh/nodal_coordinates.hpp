
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include "io/json_forward.hpp"

namespace neon
{
/// nodal_coordinates stores the list of coordinates of a discretized geometry
class nodal_coordinates
{
public:
    nodal_coordinates() = default;

    /// Construct with a list of coordinates
    explicit nodal_coordinates(matrix3x const coordinates);

    /// Construct with a list of coordinates in json format
    explicit nodal_coordinates(json const& mesh_file);

    /// \return number of nodal points
    [[nodiscard]] auto size() const noexcept { return X.cols(); }

    [[nodiscard]] matrix3x const& coordinates() const noexcept { return X; }

    /// \return the coordinates using fancy indexing
    template <typename indices_type>
    [[nodiscard]] auto coordinates(indices_type const local_node_view) const noexcept
    {
        return X(Eigen::all, local_node_view);
    }

protected:
    /// Reference configuration
    matrix3x X;
};
}
