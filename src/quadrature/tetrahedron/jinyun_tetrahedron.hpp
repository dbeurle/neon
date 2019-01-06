
#pragma once

#include "quadrature/tetrahedron/tetrahedron_quadrature.hpp"

namespace neon
{
class jinyun_tetrahedron : public tetrahedron_quadrature
{
public:
    explicit jinyun_tetrahedron(int const minimum_degree);
};
}
