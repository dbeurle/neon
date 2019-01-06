
#pragma once

#include "quadrature/numerical_quadrature.hpp"
#include "mesh/element_topology.hpp"
#include "io/json_forward.hpp"

namespace neon
{
std::unique_ptr<line_quadrature> make_line_quadrature(element_topology const topology,
                                                      int const minimum_degree,
                                                      json const& quadrature_options);

std::unique_ptr<surface_quadrature> make_surface_quadrature(element_topology const topology,
                                                            int const minimum_degree,
                                                            json const& quadrature_options);

std::unique_ptr<volume_quadrature> make_volume_quadrature(element_topology const topology,
                                                          int const minimum_degree,
                                                          json const& quadrature_options);

template <typename Quadrature>
class quadrature_factory;

/// Template specialisation for creating a volume quadrature
template <>
class quadrature_factory<volume_quadrature>
{
public:
    static std::unique_ptr<volume_quadrature> make(element_topology const topology,
                                                   int const minimum_degree,
                                                   json const& quadrature_options)
    {
        return make_volume_quadrature(topology, minimum_degree, quadrature_options);
    }
};

/// Template specialisation for creating a surface quadrature
template <>
class quadrature_factory<surface_quadrature>
{
public:
    static std::unique_ptr<surface_quadrature> make(element_topology const topology,
                                                    int const minimum_degree,
                                                    json const& quadrature_options)
    {
        return make_surface_quadrature(topology, minimum_degree, quadrature_options);
    }
};

/// Template specialisation for creating a line quadrature
template <>
class quadrature_factory<line_quadrature>
{
public:
    static std::unique_ptr<line_quadrature> make(element_topology const topology,
                                                 int const minimum_degree,
                                                 json const& quadrature_options)
    {
        return make_line_quadrature(topology, minimum_degree, quadrature_options);
    }
};

}
