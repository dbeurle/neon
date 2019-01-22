
#define EIGEN_DONT_PARALLELIZE

#include "latin_submesh.hpp"

#include "mesh/material_coordinates.hpp"
#include "numeric/mechanics"
#include "mesh/projection/recovery.hpp"

#include <range/v3/view/transform.hpp>
#include <tbb/parallel_for.h>

namespace neon::mechanics::solid
{
latin_submesh::~latin_submesh() = default;

latin_submesh::latin_submesh(latin_submesh&&) = default;

latin_submesh& latin_submesh::operator=(latin_submesh&&) = default;

auto latin_submesh::incremental_latin_internal_force(std::int32_t const element,
                                                     double const latin_search_direction) const
    -> vector const&
{
    auto const& x = coordinates->current_configuration(local_node_view(element));

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

    bilinear_gradient
        .integrate(Eigen::Map<row_matrix>(f_int_latin.data(), nodes_per_element(), dofs_per_node()),
                   [&](auto const& N_dN, auto const index) -> matrix {
                       auto const& [N, dN] = N_dN;

                       matrix3 const jacobian = local_deformation_gradient(dN, x);

                       matrix3 const& cauchy_stress_local = cauchy_stresses_local[view(element, index)];

                       matrix3 const&
                           cauchy_stress_old_local = cauchy_stresses_old_local[view(element, index)];

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

    return f_int_latin;
}
}
