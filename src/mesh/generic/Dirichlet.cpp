
#include "Dirichlet.hpp"

namespace neon
{
Dirichlet::Dirichlet(local_indices dofs, json const& times, json const& loads)
    : Boundary(times, loads), dofs(dofs)
{
}
}
