
#pragma once

#include "mesh/nodal_coordinates.hpp"

namespace neon
{
class material_coordinates : public nodal_coordinates
{
public:
    /// Construct this class using a set of initial coordinates
    material_coordinates(matrix3x const& initial_coordinates);

    /// \return element reference configuration based on the local node numbers
    template <typename indices_type>
    auto const initial_configuration(indices_type const local_nodes) const
    {
        return X(Eigen::all, local_nodes);
    }

    /// \return element current configuration based on the local node numbers
    template <typename indices_type>
    auto const current_configuration(indices_type const local_nodes) const
    {
        return x(Eigen::all, local_nodes);
    }

    /// \param u - displacement vector from initial configuration (x,y,z...)
    void update_current_configuration(vector const& u);

    /// \param u - displacement vector from initial configuration (x,y...)
    void update_current_xy_configuration(vector const& u);

    [[nodiscard]] vector displacement() const;

protected:
    /// Current configuration
    matrix3x x;
};
}
