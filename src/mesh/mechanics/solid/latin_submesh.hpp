
#pragma once

/// @file

#define EIGEN_DONT_PARALLELIZE

#include "exceptions.hpp"

#include "mesh/material_coordinates.hpp"
#include "numeric/gradient_operator.hpp"
#include "numeric/mechanics"

#include <termcolor/termcolor.hpp>

#include <tbb/parallel_for.h>

#include <cfenv>
#include <chrono>

#include <range/v3/view/transform.hpp>

#include "mesh/mechanics/solid/submesh.hpp"

namespace neon
{
namespace mechanics::solid
{
/// latin_submesh provides the element local routines for computing the system
/// components for a three-dimensional continuum mechanics discretisation.
class latin_submesh : public mechanics::solid::submesh
{
public:
    // use the base class constructor
    using mechanics::solid::submesh::submesh;

    /**
     * Compute the incremental latin internal force vector, infinite/vertical search direction,
     using the formula
     * \f{align*}{
     * f_{i} &= \int_{V} B^{T} (\sigma_{i}-\hat{\sigma}_{i}) dV
     * \f}
     * \return internal element force
     */
    [[nodiscard]] std::pair<index_view, vector const&> incremental_latin_internal_force(
        std::int32_t const element,
        double const latin_search_direction) const;
};

inline std::pair<index_view, vector const&> latin_submesh::incremental_latin_internal_force(
    std::int32_t const element,
    double const latin_search_direction) const
{
    matrix3x const& x = coordinates->current_configuration(local_node_view(element));

    auto const& cauchy_stresses_local = variables->get(variable::second::cauchy_stress);

    auto const& cauchy_stresses_old_local = variables->get_old(variable::second::cauchy_stress);

    auto const& tangent_operators = variables->get(variable::fourth::tangent_operator);

    // Compute the linear strain gradient from the displacement gradient
    auto strains_old_global = variables->get(variable::second::displacement_gradient)
                              | ranges::view::transform(
                                    [](auto const& H) { return 0.5 * (H + H.transpose()); });
    auto strains_old_local = variables->get_old(variable::second::displacement_gradient)
                             | ranges::view::transform(
                                   [](auto const& H) { return 0.5 * (H + H.transpose()); });

    static vector f_int_latin(nodes_per_element() * dofs_per_node());
    f_int_latin.setZero();

    sf->quadrature()
        .integrate_inplace(Eigen::Map<row_matrix>(f_int_latin.data(),
                                                  nodes_per_element(),
                                                  dofs_per_node()),
                           [&](auto const& N_dN, auto const index) -> matrix {
                               auto const& [N, dN] = N_dN;

                               matrix3 const jacobian = local_deformation_gradient(dN, x);

                               matrix3 const&
                                   cauchy_stress_local = cauchy_stresses_local[view(element, index)];

                               matrix3 const&
                                   cauchy_stress_old_local = cauchy_stresses_old_local[view(element,
                                                                                            index)];

                               matrix3 const& strain_old_global = strains_old_global[view(element, index)];
                               matrix3 const& strain_old_local = strains_old_local[view(element, index)];

                               matrix3 const&
                                   cauchy_stress_global = cauchy_stress_old_local
                                                          + latin_search_direction
                                                                * voigt::kinetic::from(
                                                                      tangent_operators[view(element, index)]
                                                                      * voigt::kinematic::to(
                                                                            strain_old_global
                                                                            - strain_old_local));

                               // symmetric gradient operator
                               matrix const Bt = dN * jacobian.inverse();

                               // return minus_residual = latin_residual following the LATIN implicit scheme
                               return Bt * (cauchy_stress_global - cauchy_stress_local)
                                      * jacobian.determinant();
                           });

    return {local_dof_view(element), f_int_latin};
}
}
}
