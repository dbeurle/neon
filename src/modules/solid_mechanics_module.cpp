
#include "modules/solid_mechanics_module.hpp"

#include "solver/eigen/eigen_solver_factory.hpp"

namespace neon
{
template <typename matrix_type>
solid_mechanics_module<matrix_type>::solid_mechanics_module(basic_mesh const& mesh,
                                                            json const& material,
                                                            json const& simulation)
    : fem_mesh(mesh,
               material,
               simulation["meshes"].front(),
               simulation["time"]["increments"]["initial"]),
      fem_matrix(fem_mesh, simulation)
{
}

template class solid_mechanics_module<mechanics::static_matrix<mechanics::solid::mesh>>;
template class solid_mechanics_module<mechanics::latin_matrix<mechanics::solid::mesh>>;

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
                           simulation["time"]["increments"]["initial"]},
                 make_eigen_solver(simulation["eigen_solver"])}
{
}
}
