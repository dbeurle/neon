
#include "SurfaceBoundary.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/Submesh.hpp"

#include <json/value.h>

namespace neon::diffusion
{
boundary_mesh::boundary_mesh(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                             std::vector<Submesh> const& submeshes,
                             Json::Value const& boundary,
                             Json::Value const& mesh_data)
{
    if (auto const& type = boundary["Type"].asString(); type == "HeatFlux")
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
            Json::Value heat_flux;
            for (auto i = 0; i < boundary["HeatTransferCoefficient"].size(); ++i)
            {
                heat_flux.append(Json::Value{boundary["HeatTransferCoefficient"][i].asDouble()
                                             * boundary["AmbientTemperature"][i].asDouble()});
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
