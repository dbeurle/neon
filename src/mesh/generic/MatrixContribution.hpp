
#pragma once

#include "Boundary.hpp"

namespace neon
{
/**
 * MatrixContribution is a base class for including boundary conditions which
 * contribute to the stiffness matrix in the finite element discretisation.
 */
class MatrixContribution : public Boundary
{
public:
    using Boundary::Boundary;

    /** @return an element external force vector for a given element */
    [[nodiscard]] virtual std::tuple<local_indices const&, matrix> external_stiffness(
        int const element, double const load_factor) const = 0;
};
}
