
#include "SurfaceBoundary.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/Submesh.hpp"

#include "io/json.hpp"

namespace neon::diffusion
{
boundary_mesh::boundary_mesh(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                             std::vector<Submesh> const& submeshes,
                             json const& boundary,
                             json const& mesh_data)
{
    if (auto const& type = boundary["Type"].get<std::string>(); type == "HeatFlux")
    {
        for (auto const& mesh : submeshes)
        {
            load_boundaries.emplace_back(make_surface_interpolation(mesh.topology(), mesh_data),
                                         mesh.connectivities(),
                                         mesh.connectivities(),
                                         material_coordinates,
                                         boundary["Time"],
                                         boundary["Value"]);
        }
    }
    else if (type == "NewtonCooling")
    {
        for (auto const& mesh : submeshes)
        {
            // Create the heat flux from the heat transfer coefficient and the
            // ambient temperature
            json heat_flux;
            for (auto i = 0; i < boundary["HeatTransferCoefficient"].size(); ++i)
            {
                heat_flux.emplace_back(boundary["HeatTransferCoefficient"][i].get<double>()
                                       * boundary["AmbientTemperature"][i].get<double>());
            }

            stiffness_load_boundaries.emplace_back(make_surface_interpolation(mesh.topology(),
                                                                              mesh_data),
                                                   mesh.connectivities(),
                                                   mesh.connectivities(),
                                                   material_coordinates,
                                                   boundary["Time"],
                                                   heat_flux,
                                                   boundary["HeatTransferCoefficient"]);
        }
    }
}
}
