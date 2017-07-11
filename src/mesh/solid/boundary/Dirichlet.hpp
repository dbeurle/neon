
#pragma once

#include "numeric/DenseTypes.hpp"

#include <tuple>

namespace neon::solid
{
class Dirichlet
{
public:
    Dirichlet(List dofs, double const value) : dof_value(std::make_tuple(dofs, value)) {}

    auto const& view() const { return dof_value; }

    auto const& dof_view() const { return std::get<0>(dof_value); }

    auto const& value_view() const { return std::get<1>(dof_value); }

    void update_value(double const value) { std::get<1>(dof_value) = value; }

protected:
    std::tuple<List, double> dof_value;
};
}
