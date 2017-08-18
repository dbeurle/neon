
#include "InterpolationFactory.hpp"

#include "interpolations/Hexahedron8.hpp"
#include "interpolations/Quadrilateral4.hpp"
#include "interpolations/Quadrilateral8.hpp"
#include "interpolations/Tetrahedron10.hpp"
#include "interpolations/Tetrahedron4.hpp"
#include "interpolations/Triangle3.hpp"
#include "interpolations/Triangle6.hpp"

#include <stdexcept>

#include <json/value.h>

namespace neon::solid
{
bool is_reduced_integration(Json::Value const& simulation_data)
{
    return simulation_data["ElementOptions"]["Quadrature"].empty()
               ? false
               : simulation_data["ElementOptions"]["Quadrature"].asString() == "Reduced";
}

/** Factory method for the three dimensional shape functions */
std::unique_ptr<VolumeInterpolation> make_volume_interpolation(ElementTopology const topology,
                                                               Json::Value const& simulation_data)
{
    if (!simulation_data.isMember("ElementOptions"))
    {
        throw std::runtime_error("Missing \"Part\": \"ElementOptions\"");
    }

    auto is_reduced = is_reduced_integration(simulation_data);

    switch (topology)
    {
        case ElementTopology::Hexahedron8:
        {
            return std::make_unique<Hexahedron8>(is_reduced ? HexahedronQuadrature::Rule::OnePoint
                                                            : HexahedronQuadrature::Rule::EightPoint);
        }
        case ElementTopology::Tetrahedron4:
        {
            return std::make_unique<Tetrahedron4>(TetrahedronQuadrature::Rule::OnePoint);
        }
        case ElementTopology::Tetrahedron10:
        {
            return std::make_unique<Tetrahedron10>(is_reduced
                                                       ? TetrahedronQuadrature::Rule::OnePoint
                                                       : TetrahedronQuadrature::Rule::FourPoint);
        }
        case ElementTopology::Prism6:
        case ElementTopology::Pyramid5:
        case ElementTopology::Tetrahedron20:
        case ElementTopology::Tetrahedron35:
        case ElementTopology::Tetrahedron56:
        case ElementTopology::Hexahedron64:
        case ElementTopology::Hexahedron125:
        case ElementTopology::Hexahedron27:
        case ElementTopology::Prism18:
        case ElementTopology::Pyramid14:
        case ElementTopology::Hexahedron20:
        case ElementTopology::Prism15:
        case ElementTopology::Pyramid13:
        default:
            throw std::runtime_error("Element shape not implemented for continuum "
                                     "simulations\n");
            break;
    }
    return nullptr;
}

std::unique_ptr<SurfaceInterpolation> make_surface_interpolation(ElementTopology const topology,
                                                                 Json::Value const& simulation_data)
{
    if (!simulation_data.isMember("ElementOptions"))
    {
        throw std::runtime_error("Missing \"Part\": \"ElementOptions\"");
    }

    auto is_reduced = is_reduced_integration(simulation_data);

    switch (topology)
    {
        case ElementTopology::Quadrilateral4:
        {
            return std::make_unique<Quadrilateral4>(is_reduced
                                                        ? QuadrilateralQuadrature::Rule::OnePoint
                                                        : QuadrilateralQuadrature::Rule::FourPoint);
        }
        case ElementTopology::Triangle3:
        {
            return std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint);
        }
        case ElementTopology::Triangle6:
        {
            return std::make_unique<Triangle6>(is_reduced ? TriangleQuadrature::Rule::OnePoint
                                                          : TriangleQuadrature::Rule::ThreePoint);
        }
        case ElementTopology::Quadrilateral8:
        {
            return std::make_unique<Quadrilateral8>(is_reduced
                                                        ? QuadrilateralQuadrature::Rule::FourPoint
                                                        : QuadrilateralQuadrature::Rule::NinePoint);
        }
        default:
            throw std::runtime_error("Surface element shape not implemented for "
                                     "continuum simulations\n");
            break;
    }
    return nullptr;
}
}
