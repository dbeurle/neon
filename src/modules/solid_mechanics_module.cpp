
#include "modules/solid_mechanics_module.hpp"

namespace neon
{
solid_mechanics_module::solid_mechanics_module(basic_mesh const& mesh,
                                           json const& material,
                                           json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0]), fem_matrix(fem_mesh, simulation)
{
}
}
