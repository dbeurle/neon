
#include "surface_boundary.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/basic_submesh.hpp"

#include "io/json.hpp"

namespace neon::diffusion::reaction
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
}
}
