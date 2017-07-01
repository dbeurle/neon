
#pragma once

namespace neon::solid
{
class Dirichlet
{
public:
    Dirichlet(List dofs, double const value) : dofs(dofs), value(value) {}

    auto const& prescribed_dofs() const { return dofs; }
    auto prescribed_value() const { return value; }

protected:
    List dofs;
    double value;
};
}
