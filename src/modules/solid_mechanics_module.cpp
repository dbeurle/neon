
#include "modules/solid_mechanics_module.hpp"

namespace neon
{
solid_mechanics_module::solid_mechanics_module(basic_mesh const& mesh,
                                               json const& material,
                                               json const& simulation)
    : fem_mesh(mesh,
               material,
               simulation["meshes"].front(),
               simulation["time"]["increments"]["initial"]),
      fem_matrix(fem_mesh, simulation)
{
}

linear_buckling_module::linear_buckling_module(basic_mesh const& mesh,
                                               json const& material,
                                               json const& simulation)
    : fem_mesh(mesh,
               material,
               simulation["meshes"].front(),
               simulation["time"]["increments"]["initial"]),
      fem_matrix(fem_mesh)
{
}

natural_frequency_module::natural_frequency_module(basic_mesh const& mesh,
                                                   json const& material,
                                                   json const& simulation)
    : fem_matrix{mesh_type{mesh,
                           material,
                           simulation["meshes"].front(),
                           simulation["time"]["increments"]["initial"]}}
{
}
}
