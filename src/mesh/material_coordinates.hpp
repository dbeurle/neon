
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
        return X(Eigen::placeholders::all, local_nodes);
    }

    /// \return element current configuration based on the local node numbers
    template <typename indices_type>
    auto const current_configuration(indices_type const local_nodes) const
    {
        return x(Eigen::placeholders::all, local_nodes);
    }

    /// \param u - displacement vector from initial configuration (x,y,z...)
    void update_current_configuration(vector const& u);

    /// \param u - displacement vector from initial configuration (x,y...)
    void update_current_xy_configuration(vector const& u);

    [[nodiscard]] matrix3x displacement() const;

    [[nodiscard]] vector displacement_vector() const;

protected:
    /// Current configuration
    matrix3x x;
};

// template <typename traits>
// class material_coordinates_test : public mesh_coordinates<traits>
// {
//     using coordinate_t = typename mesh_coordinates<traits>::coordinate_t;
//
//     using base_type = mesh_coordinates<traits>;
//
//     using base_type::fixed_size;
//
// public:
//     /** \return element reference configuration based on the local node numbers*/
//     auto const initial_configuration(index_view local_nodes) const
//     {
//         return X(Eigen::placeholders::all, local_nodes);
//     }
//
//     /** @return element current configuration based on the local node numbers*/
//     auto const current_configuration(index_view local_nodes) const
//     {
//         return x(Eigen::placeholders::all, local_nodes);
//     }
//
//     /** @param u - displacement vector from initial configuration */
//     void update_current_configuration(vector const& u);
//
//     [[nodiscard]] coordinate_t displacement() const { return x - X; }
//
//     protected : using mesh_coordinates<traits>::X; //!< Initial configuration
//     coordinate_t x;                                //!< Current configuration
// };
//
// template <typename traits>
// void material_coordinates<traits>::update_current_configuration(vector const& u)
// {
//     for (auto i = 0; i < fixed_size; i++)
//     {
//         x.row(i) = X.row(i) + u.transpose()(Eigen::seq(i, u.size() - 1, fixed_size));
//     }
// }
}
