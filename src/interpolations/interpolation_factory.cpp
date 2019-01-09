
#include "interpolation_factory.hpp"

#include "interpolations/shape_function.hpp"

#include "interpolations/hexahedron.hpp"
#include "interpolations/prism.hpp"
#include "interpolations/tetrahedron.hpp"
#include "interpolations/pyramid.hpp"

#include "interpolations/quadrilateral.hpp"
#include "interpolations/triangle.hpp"

#include "interpolations/line.hpp"

#include <stdexcept>

namespace neon
{
std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology)
{
    switch (topology)
    {
        case element_topology::hexahedron8:
        {
            return std::make_unique<hexahedron8>();
        }
        case element_topology::hexahedron20:
        {
            return std::make_unique<hexahedron20>();
        }
        case element_topology::hexahedron27:
        {
            return std::make_unique<hexahedron27>();
        }
        case element_topology::tetrahedron4:
        {
            return std::make_unique<tetrahedron4>();
        }
        case element_topology::tetrahedron10:
        {
            return std::make_unique<tetrahedron10>();
        }
        case element_topology::prism6:
        {
            return std::make_unique<prism6>();
        }
        case element_topology::prism15:
        {
            return std::make_unique<prism15>();
        }
        case element_topology::pyramid5:
        {
            return std::make_unique<pyramid5>();
        }
        case element_topology::pyramid13:
        {
            return std::make_unique<pyramid13>();
        }
        case element_topology::pyramid14:
        case element_topology::prism18:
        case element_topology::tetrahedron20:
        case element_topology::tetrahedron35:
        case element_topology::tetrahedron56:
        case element_topology::hexahedron64:
        case element_topology::hexahedron125:
        default:
        {
            throw std::runtime_error("Volume interpolation number "
                                     + std::to_string(static_cast<int>(topology))
                                     + " not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology)
{
    switch (topology)
    {
        case element_topology::quadrilateral4:
        {
            return std::make_unique<quadrilateral4>();
        }
        case element_topology::triangle3:
        {
            return std::make_unique<triangle3>();
        }
        case element_topology::triangle6:
        {
            return std::make_unique<triangle6>();
        }
        case element_topology::quadrilateral8:
        {
            return std::make_unique<quadrilateral8>();
        }
        case element_topology::quadrilateral9:
        {
            return std::make_unique<quadrilateral9>();
        }
        default:
        {
            throw std::runtime_error("Surface element shape not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology)
{
    switch (topology)
    {
        case element_topology::line2:
        {
            return std::make_unique<line2>();
        }
        case element_topology::line3:
        {
            return std::make_unique<line3>();
        }
        default:
            throw std::runtime_error("Line element shape not implemented");
            break;
    }
    return nullptr;
}
}
