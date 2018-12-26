
#define EIGEN_DONT_PARALLELIZE

#include "latin_submesh.hpp"

#include "exceptions.hpp"

#include "constitutive/constitutive_model_factory.hpp"
#include "interpolations/interpolation_factory.hpp"
#include "material/material_property.hpp"
#include "mesh/material_coordinates.hpp"
#include "numeric/gradient_operator.hpp"
#include "numeric/mechanics"
#include "mesh/dof_allocator.hpp"
#include "traits/mechanics.hpp"

#include <termcolor/termcolor.hpp>

#include <tbb/parallel_for.h>

#include <cfenv>
#include <chrono>

namespace neon::mechanics::solid
{
std::pair<index_view, vector const&> latin_submesh::incremental_latin_internal_force(
    std::int32_t const element) const
{
    auto const& x = coordinates->current_configuration(local_node_view(element));

    auto const& cauchy_stresses_local = variables->get(variable::second::cauchy_stress);

    auto const& cauchy_stresses_old_global = variables->get_old(variable::second::cauchy_stress);

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

                               matrix3 const& cauchy_stress_old_global = cauchy_stresses_old_global
                                   [view(element, index)];

                               // symmetric gradient operator
                               matrix const Bt = dN * jacobian.inverse();

                               return Bt * (cauchy_stress_old_global - cauchy_stress_local)
                                      * jacobian.determinant();
                           });

    return {local_dof_view(element), f_int_latin};
}
}
