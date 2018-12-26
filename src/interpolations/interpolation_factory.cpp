
#include "interpolation_factory.hpp"

#include "interpolations/shape_function.hpp"

#include "interpolations/hexahedron.hpp"
#include "interpolations/prism.hpp"
#include "interpolations/tetrahedron.hpp"

#include "interpolations/quadrilateral.hpp"
#include "interpolations/triangle.hpp"

#include "interpolations/line.hpp"

#include <stdexcept>

#include <io/json.hpp>

namespace neon
{
static bool is_reduced_integration(json const& simulation_data)
{
    return simulation_data["element_options"]["quadrature"].is_null()
               ? false
               : simulation_data["element_options"]["quadrature"].get<std::string>() == "reduced";
}

static void check_element_options(json const& simulation_data)
{
    if (simulation_data.find("element_options") == end(simulation_data))
    {
        throw std::runtime_error("Missing \"part\": {\"element_options\" : {}}");
    }
}

std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology,
                                                                json const& simulation_data)
{
    check_element_options(simulation_data);

    auto const is_reduced = is_reduced_integration(simulation_data);

    switch (topology)
    {
        case element_topology::hexahedron8:
        {
            return std::make_unique<hexahedron8>(is_reduced ? hexahedron_quadrature::point::one
                                                            : hexahedron_quadrature::point::eight);
        }
        case element_topology::hexahedron20:
        {
            return std::make_unique<hexahedron20>(is_reduced
                                                      ? hexahedron_quadrature::point::eight
                                                      : hexahedron_quadrature::point::twentyseven);
        }
        case element_topology::hexahedron27:
        {
            return std::make_unique<hexahedron27>(is_reduced
                                                      ? hexahedron_quadrature::point::eight
                                                      : hexahedron_quadrature::point::twentyseven);
        }
        case element_topology::tetrahedron4:
        {
            return std::make_unique<tetrahedron4>(tetrahedron_quadrature::point::one);
        }
        case element_topology::tetrahedron10:
        {
            return std::make_unique<tetrahedron10>(is_reduced ? tetrahedron_quadrature::point::one
                                                              : tetrahedron_quadrature::point::four);
        }
        case element_topology::prism6:
        {
            return std::make_unique<prism6>(is_reduced ? prism_quadrature::point::one
                                                       : prism_quadrature::point::six);
        }
        case element_topology::prism15:
        {
            return std::make_unique<prism15>(is_reduced ? prism_quadrature::point::six
                                                        : prism_quadrature::point::nine);
        }
        case element_topology::prism18:
        case element_topology::pyramid5:
        case element_topology::pyramid13:
        case element_topology::pyramid14:
        case element_topology::tetrahedron20:
        case element_topology::tetrahedron35:
        case element_topology::tetrahedron56:
        case element_topology::hexahedron64:
        case element_topology::hexahedron125:
        default:
            throw std::runtime_error("Volume interpolation "
                                     + std::to_string(static_cast<int>(topology))
                                     + " not implemented");
            break;
    }
    return nullptr;
}

std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology,
                                                                  json const& simulation_data)
{
    check_element_options(simulation_data);

    auto is_reduced = is_reduced_integration(simulation_data);

    switch (topology)
    {
        case element_topology::quadrilateral4:
        {
            return std::make_unique<quadrilateral4>(is_reduced
                                                        ? quadrilateral_quadrature::point::one
                                                        : quadrilateral_quadrature::point::four);
        }
        case element_topology::triangle3:
        {
            return std::make_unique<triangle3>(triangle_quadrature::point::one);
        }
        case element_topology::triangle6:
        {
            return std::make_unique<triangle6>(is_reduced ? triangle_quadrature::point::one
                                                          : triangle_quadrature::point::three);
        }
        case element_topology::quadrilateral8:
        {
            return std::make_unique<quadrilateral8>(is_reduced
                                                        ? quadrilateral_quadrature::point::four
                                                        : quadrilateral_quadrature::point::nine);
        }
        case element_topology::quadrilateral9:
        {
            return std::make_unique<quadrilateral9>(is_reduced
                                                        ? quadrilateral_quadrature::point::four
                                                        : quadrilateral_quadrature::point::nine);
        }
        default:
        {
            throw std::runtime_error("Surface element shape not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology,
                                                            json const& simulation_data)
{
    check_element_options(simulation_data);

    auto is_reduced = is_reduced_integration(simulation_data);

    switch (topology)
    {
        case element_topology::line2:
        {
            return std::make_unique<line2>(line_quadrature::point::one);
        }
        case element_topology::line3:
        {
            return std::make_unique<line3>(is_reduced ? line_quadrature::point::one
                                                      : line_quadrature::point::two);
        }
        default:
            throw std::runtime_error("Line element shape not implemented");
            break;
    }
    return nullptr;
}
}
