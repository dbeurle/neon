
#include "quadrature_factory.hpp"

#include "quadrature/hexahedron/hexahedron_quadrature.hpp"

#include "quadrature/line/gauss_legendre_line.hpp"

#include "quadrature/prism/felippa_prism.hpp"
#include "quadrature/prism/kubatko_prism.hpp"

#include "quadrature/pyramid/bedrosian_pyramid.hpp"

#include "quadrature/tetrahedron/jinyun_tetrahedron.hpp"
#include "quadrature/tetrahedron/witherden_vincent_tetrahedron.hpp"

#include "quadrature/quadrilateral/quadrilateral_quadrature.hpp"

#include "quadrature/triangle/cowper_triangle.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon
{
std::unique_ptr<volume_quadrature> make_volume_quadrature(element_topology const topology,
                                                          int const minimum_degree,
                                                          json const& quadrature_options)
{
    switch (topology)
    {
        case element_topology::hexahedron8:
        case element_topology::hexahedron20:
        case element_topology::hexahedron27:
        case element_topology::hexahedron64:
        case element_topology::hexahedron125:
        {
            if (quadrature_options.find("hexahedron") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<hexahedron_quadrature>(minimum_degree);
        }
        case element_topology::tetrahedron4:
        case element_topology::tetrahedron10:
        case element_topology::tetrahedron20:
        case element_topology::tetrahedron35:
        case element_topology::tetrahedron56:
        {
            if (quadrature_options.find("tetrahedron") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<witherden_vincent_tetrahedron>(minimum_degree);
        }
        case element_topology::prism6:
        case element_topology::prism15:
        case element_topology::prism18:
        {
            if (quadrature_options.find("prism") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<felippa_prism>(minimum_degree);
        }
        case element_topology::pyramid5:
        case element_topology::pyramid13:
        case element_topology::pyramid14:
        {
            if (quadrature_options.find("pyramid") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<bedrosian_pyramid>(minimum_degree);
        }
        default:
        {
            throw std::runtime_error("Volume quadrature for element number "
                                     + std::to_string(static_cast<int>(topology))
                                     + " is not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<surface_quadrature> make_surface_quadrature(element_topology const topology,
                                                            int const minimum_degree,
                                                            json const& quadrature_options)
{
    switch (topology)
    {
        case element_topology::quadrilateral4:
        case element_topology::quadrilateral8:
        case element_topology::quadrilateral9:
        {
            if (quadrature_options.find("quadrilateral") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<quadrilateral_quadrature>(minimum_degree);
        }
        case element_topology::triangle3:
        case element_topology::triangle6:
        {
            if (quadrature_options.find("triangle") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<cowper_triangle>(minimum_degree);
        }
        default:
        {
            throw std::runtime_error("Surface element shape not implemented");
            break;
        }
    }
    return nullptr;
}

std::unique_ptr<line_quadrature> make_line_quadrature(element_topology const topology,
                                                      int const minimum_degree,
                                                      json const& quadrature_options)
{
    switch (topology)
    {
        case element_topology::line2:
        case element_topology::line3:
        {
            if (quadrature_options.find("line") == end(quadrature_options))
            {
                // provide a default
            }
            return std::make_unique<gauss_legendre_line>(minimum_degree);
        }
        default:
            throw std::runtime_error("Line element shape not implemented");
            break;
    }
    return nullptr;
}
}
