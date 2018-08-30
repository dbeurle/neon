
#include "modules/beam_module.hpp"

namespace neon
{
beam_module::beam_module(basic_mesh const& mesh,
                         json const& material_data,
                         json const& simulation_data,
                         std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store)
    : fem_mesh(mesh,
               material_data,
               simulation_data["meshes"].front(),
               simulation_data["time"]["increments"]["initial"],
               profile_store),
      fem_matrix(fem_mesh, simulation_data)
{
}
}
