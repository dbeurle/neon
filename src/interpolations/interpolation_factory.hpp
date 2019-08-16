
#pragma once

/// @file

#include "shape_function_forward.hpp"
#include "mesh/element_topology.hpp"

#include <memory>

namespace neon
{
/// Factory method for the three dimensional shape functions
/// \param topology Element topology
std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology);

template <typename Interpolation>
class interpolation_factory;

/// Template specialisation for creating a volume interpolation
template <>
class interpolation_factory<volume_interpolation>
{
public:
    static std::unique_ptr<volume_interpolation> make(element_topology const topology)
    {
        return make_volume_interpolation(topology);
    }
};

/// Template specialisation for creating a surface interpolation
template <>
class interpolation_factory<surface_interpolation>
{
public:
    static std::unique_ptr<surface_interpolation> make(element_topology const topology)
    {
        return make_surface_interpolation(topology);
    }
};

/// Template specialisation for creating a line interpolation
template <>
class interpolation_factory<line_interpolation>
{
public:
    static std::unique_ptr<line_interpolation> make(element_topology const topology)
    {
        return make_line_interpolation(topology);
    }
};

}
