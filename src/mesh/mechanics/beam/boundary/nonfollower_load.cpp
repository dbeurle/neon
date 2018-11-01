
#include "nonfollower_load.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "io/json.hpp"
#include "traits/mechanics.hpp"

namespace neon::mechanics::beam
{
nonfollower_load_boundary::nonfollower_load_boundary(
    std::shared_ptr<material_coordinates>& material_coordinates,
    std::vector<basic_submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary_data,
    std::unordered_map<std::string, int> const& dof_table,
    double const generate_time_step)
{
    if (std::string const& type = boundary_data["type"]; type == "moment" || type == "force")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.find(dof_name) != boundary_data.end())
            {
                for (auto const& mesh : submeshes)
                {
                    // create linear dof indices
                    auto node_indices = mesh.unique_node_indices();

                    // Offset the degrees of freedom on the boundary
                    std::transform(begin(node_indices),
                                   end(node_indices),
                                   begin(node_indices),
                                   [&, dof_offset = std::ref(dof_offset)](auto const node) {
                                       return node * 6 + dof_offset + (type == "force" ? 0 : 3);
                                   });

                    nodal_values.emplace_back(node_indices, boundary_data, dof_name, generate_time_step);
                }
            }
        }
    }
    else
    {
        throw std::domain_error("Specify a valid boundary type \"force\" or \"moment\"");
    }
}
}
