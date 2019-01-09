
#pragma once

#include "quadrature/tetrahedron/tetrahedron_quadrature.hpp"

namespace neon
{
class witherden_vincent_tetrahedron : public tetrahedron_quadrature
{
public:
    explicit witherden_vincent_tetrahedron(int const minimum_degree);
};
}
