
#include "Dirichlet.hpp"

namespace neon
{
Dirichlet::Dirichlet(List dofs, json const& times, json const& loads)
    : Boundary(times, loads), dofs(dofs)
{
}
}
