
#include "Dirichlet.hpp"

namespace neon::solid
{
Dirichlet::Dirichlet(List dofs, double const value, bool const is_ramped)
    : dofs(dofs), is_ramped(is_ramped), value_new(value)
{
}

void Dirichlet::update_value(double const value)
{
    value_old = value_new;
    value_new = value;
}

void Dirichlet::inherit_from_last()
{
    value_old = value_new;
    is_ramped = false;
}
}
