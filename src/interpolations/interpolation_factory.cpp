
#include "interpolation_factory.hpp"

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
bool is_reduced_integration(json const& mesh_data)
{
    return mesh_data["ElementOptions"]["Quadrature"].is_null()
               ? false
               : mesh_data["ElementOptions"]["Quadrature"].get<std::string>() == "Reduced";
}

void check_element_options(json const& mesh_data)
{
    if (!mesh_data.count("ElementOptions"))
    {
        throw std::runtime_error("Missing \"Part\": \"ElementOptions\"");
    }
}

std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology,
                                                                json const& mesh_data)
{
    check_element_options(mesh_data);

    auto const is_reduced = is_reduced_integration(mesh_data);

    switch (topology)
    {
        case element_topology::hexahedron8:
        {
            return std::make_unique<hexahedron8>(is_reduced
                                                     ? hexahedron_quadrature::Rule::OnePoint
                                                     : hexahedron_quadrature::Rule::EightPoint);
        }
        case element_topology::hexahedron20:
        {
            return std::make_unique<hexahedron20>(
                is_reduced ? hexahedron_quadrature::Rule::EightPoint
                           : hexahedron_quadrature::Rule::TwentySevenPoint);
        }
        case element_topology::hexahedron27:
        {
            return std::make_unique<hexahedron27>(
                is_reduced ? hexahedron_quadrature::Rule::EightPoint
                           : hexahedron_quadrature::Rule::TwentySevenPoint);
        }
        case element_topology::tetrahedron4:
        {
            return std::make_unique<tetrahedron4>(tetrahedron_quadrature::Rule::OnePoint);
        }
        case element_topology::tetrahedron10:
        {
            return std::make_unique<tetrahedron10>(is_reduced
                                                       ? tetrahedron_quadrature::Rule::OnePoint
                                                       : tetrahedron_quadrature::Rule::FourPoint);
        }
        case element_topology::prism6:
        {
            return std::make_unique<prism6>(prism_quadrature::Rule::OnePoint);
        }
        case element_topology::prism15:
        {
            return std::make_unique<prism15>(is_reduced ? prism_quadrature::Rule::OnePoint
                                                        : prism_quadrature::Rule::SixPoint);
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
                                                                  json const& mesh_data)
{
    check_element_options(mesh_data);

    auto is_reduced = is_reduced_integration(mesh_data);

    switch (topology)
    {
        case element_topology::quadrilateral4:
        {
            return std::make_unique<quadrilateral4>(is_reduced
                                                        ? quadrilateral_quadrature::Rule::OnePoint
                                                        : quadrilateral_quadrature::Rule::FourPoint);
        }
        case element_topology::triangle3:
        {
            return std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint);
        }
        case element_topology::triangle6:
        {
            return std::make_unique<triangle6>(is_reduced ? triangle_quadrature::Rule::OnePoint
                                                          : triangle_quadrature::Rule::ThreePoint);
        }
        case element_topology::quadrilateral8:
        {
            return std::make_unique<quadrilateral8>(is_reduced
                                                        ? quadrilateral_quadrature::Rule::FourPoint
                                                        : quadrilateral_quadrature::Rule::NinePoint);
        }
        case element_topology::quadrilateral9:
        {
            return std::make_unique<quadrilateral9>(is_reduced
                                                        ? quadrilateral_quadrature::Rule::FourPoint
                                                        : quadrilateral_quadrature::Rule::NinePoint);
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
                                                            json const& mesh_data)
{
    check_element_options(mesh_data);

    auto is_reduced = is_reduced_integration(mesh_data);

    switch (topology)
    {
        case element_topology::line2:
        {
            return std::make_unique<line2>(line_quadrature::Rule::OnePoint);
        }
        case element_topology::line3:
        {
            return std::make_unique<line3>(is_reduced ? line_quadrature::Rule::OnePoint
                                                      : line_quadrature::Rule::TwoPoint);
        }
        default:
            throw std::runtime_error("Line element shape not implemented");
            break;
    }
    return nullptr;
}
}
