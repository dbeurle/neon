
#include "InterpolationFactory.hpp"

#include "interpolations/Hexahedron.hpp"
#include "interpolations/Quadrilateral.hpp"

#include "interpolations/Tetrahedron10.hpp"
#include "interpolations/Tetrahedron4.hpp"
#include "interpolations/Triangle3.hpp"
#include "interpolations/Triangle6.hpp"

#include "interpolations/Line.hpp"
#include "interpolations/prism.hpp"

#include <stdexcept>

#include <json/value.h>

namespace neon
{
bool is_reduced_integration(Json::Value const& mesh_data)
{
    return mesh_data["ElementOptions"]["Quadrature"].empty()
               ? false
               : mesh_data["ElementOptions"]["Quadrature"].asString() == "Reduced";
}

void check_element_options(Json::Value const& mesh_data)
{
    if (!mesh_data.isMember("ElementOptions"))
    {
        throw std::runtime_error("Missing \"Part\": \"ElementOptions\"");
    }
}

/** Factory method for the three dimensional shape functions */
std::unique_ptr<VolumeInterpolation> make_volume_interpolation(ElementTopology const topology,
                                                               Json::Value const& mesh_data)
{
    check_element_options(mesh_data);

    auto const is_reduced = is_reduced_integration(mesh_data);

    switch (topology)
    {
        case ElementTopology::Hexahedron8:
        {
            return std::make_unique<Hexahedron8>(is_reduced ? HexahedronQuadrature::Rule::OnePoint
                                                            : HexahedronQuadrature::Rule::EightPoint);
        }
        case ElementTopology::Hexahedron20:
        {
            return std::make_unique<Hexahedron20>(is_reduced
                                                      ? HexahedronQuadrature::Rule::EightPoint
                                                      : HexahedronQuadrature::Rule::TwentySevenPoint);
        }
        case ElementTopology::Hexahedron27:
        {
            return std::make_unique<Hexahedron27>(is_reduced
                                                      ? HexahedronQuadrature::Rule::EightPoint
                                                      : HexahedronQuadrature::Rule::TwentySevenPoint);
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
        {
            return std::make_unique<prism6>(PrismQuadrature::Rule::OnePoint);
        }
        case ElementTopology::Prism15:
        {
            return std::make_unique<prism15>(is_reduced ? PrismQuadrature::Rule::OnePoint
                                                        : PrismQuadrature::Rule::SixPoint);
        }
        case ElementTopology::Prism18:
        case ElementTopology::Pyramid5:
        case ElementTopology::Pyramid13:
        case ElementTopology::Pyramid14:
        case ElementTopology::Tetrahedron20:
        case ElementTopology::Tetrahedron35:
        case ElementTopology::Tetrahedron56:
        case ElementTopology::Hexahedron64:
        case ElementTopology::Hexahedron125:
        default:
            throw std::runtime_error("Volume interpolation "
                                     + std::to_string(static_cast<int>(topology))
                                     + " not implemented");
            break;
    }
    return nullptr;
}

std::unique_ptr<SurfaceInterpolation> make_surface_interpolation(ElementTopology const topology,
                                                                 Json::Value const& mesh_data)
{
    check_element_options(mesh_data);

    auto is_reduced = is_reduced_integration(mesh_data);

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
        case ElementTopology::Quadrilateral9:
        {
            return std::make_unique<Quadrilateral9>(is_reduced
                                                        ? QuadrilateralQuadrature::Rule::FourPoint
                                                        : QuadrilateralQuadrature::Rule::NinePoint);
        }
        default:
        {
            throw std::runtime_error("Surface element shape not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<LineInterpolation> make_line_interpolation(ElementTopology const topology,
                                                           Json::Value const& mesh_data)
{
    check_element_options(mesh_data);

    auto is_reduced = is_reduced_integration(mesh_data);

    switch (topology)
    {
        case ElementTopology::Line2:
        {
            return std::make_unique<Line2>(LineQuadrature::Rule::OnePoint);
        }
        case ElementTopology::Line3:
        {
            return std::make_unique<Line3>(is_reduced ? LineQuadrature::Rule::OnePoint
                                                      : LineQuadrature::Rule::TwoPoint);
        }
        default:
            throw std::runtime_error("Line element shape not implemented");
            break;
    }
    return nullptr;
}
}
