
#include "surface_boundary.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/basic_submesh.hpp"

#include "io/json.hpp"

namespace neon::diffusion
{
boundary_mesh::boundary_mesh(std::shared_ptr<material_coordinates>& material_coordinates,
                             std::vector<basic_submesh> const& submeshes,
                             json const& boundary,
                             json const& mesh_data)
{
    if (std::string const& type = boundary["type"]; type == "flux")
    {
        for (auto const& mesh : submeshes)
        {
            load_boundaries.emplace_back(make_surface_interpolation(mesh.topology(), mesh_data),
                                         mesh.all_node_indices(),
                                         mesh.all_node_indices(),
                                         material_coordinates,
                                         boundary["time"],
                                         boundary["value"]);
        }
    }
    else if (type == "newton_cooling")
    {
        for (auto const& mesh : submeshes)
        {
            // Create the heat flux from the heat transfer coefficient and the
            // ambient temperature
            json heat_flux;
            for (std::size_t i{0}; i < boundary["heat_transfer_coefficient"].size(); ++i)
            {
                heat_flux.emplace_back(boundary["heat_transfer_coefficient"][i].get<double>()
                                       * boundary["ambient_temperature"][i].get<double>());
            }

            stiffness_load_boundaries.emplace_back(make_surface_interpolation(mesh.topology(),
                                                                              mesh_data),
                                                   mesh.all_node_indices(),
                                                   mesh.all_node_indices(),
                                                   material_coordinates,
                                                   boundary["time"],
                                                   heat_flux,
                                                   boundary["heat_transfer_coefficient"]);
        }
    }
}
}
