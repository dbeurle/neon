
#include "modules/beam_module.hpp"

namespace neon
{
beam_module::beam_module(basic_mesh const& mesh, json const& material, json const& simulation_data)
    : fem_mesh(mesh,
               material,
               simulation_data["Mesh"][0],
               simulation_data["Time"]["Increments"]["Initial"]),
      fem_matrix(fem_mesh, simulation_data)
{
}
}
