
#include "modules/solid_mechanics_module.hpp"

namespace neon
{
template <typename matrix_type>
solid_mechanics_module<matrix_type>::solid_mechanics_module(basic_mesh const& mesh,
                                                            json const& material,
                                                            json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0], simulation["Time"]["Increments"]["Initial"]),
      fem_matrix(fem_mesh, simulation)
{
}
template class solid_mechanics_module<mechanics::static_matrix<mechanics::solid::mesh>>;
template class solid_mechanics_module<mechanics::latin_matrix<mechanics::solid::mesh>>;

solid_mechanics_linear_buckling_module::solid_mechanics_linear_buckling_module(basic_mesh const& mesh,
                                                                               json const& material,
                                                                               json const& simulation)
    : fem_mesh(mesh, material, simulation["Mesh"][0], simulation["Time"]["Increments"]["Initial"]),
      fem_matrix(fem_mesh)
{
}
}
