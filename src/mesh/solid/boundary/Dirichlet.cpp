
#include "Dirichlet.hpp"

namespace neon::solid
{
Dirichlet::Dirichlet(List dofs, double const prescribed_value, bool const is_load_ramped)
    : Boundary(prescribed_value, is_load_ramped), dofs(dofs)
{
}
}
