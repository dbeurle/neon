
#pragma once

#include "mesh/nodal_coordinates_v2.hpp"

namespace neon::v2
{
template <int N>
class material_coordinates : public nodal_coordinates<N>
{
public:
    using value_type = typename nodal_coordinates<N>::value_type;

public:
    explicit material_coordinates(value_type const& initial_coordinates)
        : nodal_coordinates<N>(initial_coordinates), x(initial_coordinates)
    {
    }

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

    /// \param u displacement vector from initial configuration (x,y,z...)
    void update_current_configuration(vector const& u)
    {
        for (int n = 0; n < N; n++)
        {
            x.row(n) = X.row(n) + u.transpose()(Eigen::seq(n, u.size() - 1, N));
        }
    }

    [[nodiscard]] value_type displacement() const { return x - X; }

        [[nodiscard]] vector displacement_vector() const
    {
        return Eigen::Map<vector>(displacement().data(), displacement().size());
    }

protected:
    /// Current configuration
    value_type x;

private:
    using nodal_coordinates<N>::X;
};

extern template class material_coordinates<2>;
extern template class material_coordinates<3>;
}
