
#include "femSubmesh.hpp"

#include "Exceptions.hpp"

#include "constitutive/ConstitutiveModelFactory.hpp"
#include "interpolations/InterpolationFactory.hpp"
#include "material/Material.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"
#include "numeric/Operators.hpp"

#include <cfenv>
#include <chrono>
#include <omp.h>

#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/transform.hpp>

#include <termcolor/termcolor.hpp>

namespace neon::mechanical::plane
{
femSubmesh::femSubmesh(Json::Value const& material_data,
                       Json::Value const& simulation_data,
                       std::shared_ptr<MaterialCoordinates>& material_coordinates,
                       Submesh const& submesh)
    : detail::femSubmesh<plane::femSubmesh>(submesh),
      material_coordinates(material_coordinates),
      sf(solid::make_surface_interpolation(topology(), simulation_data)),
      variables(elements() * sf->quadrature().points()),
      cm(make_constitutive_model(variables, material_data, simulation_data))
{
    // Allocate storage for the displacement gradient
    variables.add(InternalVariables::Tensor::DisplacementGradient,
                  InternalVariables::Tensor::DeformationGradient,
                  InternalVariables::Tensor::Cauchy);

    variables.add(InternalVariables::Scalar::DetF);

    // Get the old data to the undeformed configuration
    for (auto& F : variables(InternalVariables::Tensor::DeformationGradient))
    {
        F = Matrix2::Identity();
    }
    variables.commit();

    dof_list = allocate_dof_list(this->dofs_per_node(), nodal_connectivity);
}

void femSubmesh::save_internal_variables(bool const have_converged)
{
    if (have_converged)
    {
        variables.commit();
    }
    else
    {
        variables.revert();
    }
}

std::tuple<List const&, Matrix> femSubmesh::tangent_stiffness(int64 const element) const
{
    auto const x = material_coordinates->current_configuration(local_node_list(element));

    Matrix const kmat = material_tangent_stiffness(x, element);

    if (!cm->is_finite_deformation()) return {local_dof_list(element), kmat};

    Matrix const kgeo = geometric_tangent_stiffness(x, element);

    return {local_dof_list(element), kmat + kgeo};
}

std::tuple<List const&, Vector> femSubmesh::internal_force(int64 const element) const
{
    auto const x = material_coordinates->current_configuration(local_node_list(element));

    return {local_dof_list(element), internal_nodal_force(x, element)};
}

Matrix femSubmesh::geometric_tangent_stiffness(Matrix const& x, int64 const element) const
{
    // auto const& cauchy_list = variables(InternalVariables::Tensor::Cauchy);

    // auto const n = nodes_per_element();

    return Matrix::Ones(1, 1);

    // clang-format off
    // Matrix const kgeo = sf->quadrature().integrate(Matrix::Zero(n, n).eval(),
    //                                                [&](auto const& femval, auto const& l) -> Matrix {
    //                                                     auto const & [ N, rhea ] = femval;
    //
    //                                                     auto const Jacobian = local_deformation_gradient(rhea, x);
    //
    //                                                     Matrix2 const cauchy = cauchy_list[offset(element, l)];
    //
    //                                                     // Compute the symmetric gradient operator
    //                                                     auto const L = (rhea * Jacobian.inverse()).transpose();
    //
    //                                                     return L.transpose() * cauchy * L * Jacobian.determinant();
    //                                                 });
    // // clang-format on
    // return identity_expansion(kgeo, dofs_per_node());
}

Matrix femSubmesh::material_tangent_stiffness(Matrix const& x, int64 const element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    // auto const& tangent_operators = variables(InternalVariables::Matrix::TangentOperator);

    Matrix kmat = Matrix::Zero(local_dofs, local_dofs);

    return kmat;

    // return sf->quadrature().integrate(kmat, [&](auto const& femval, auto const& l) -> Matrix {
    //
    //     auto const & [N, rhea] = femval;
    //
    //     auto const& D = tangent_operators[offset(element, l)];
    //
    //     auto const Jacobian = local_deformation_gradient(rhea, x);
    //
    //     // Compute the symmetric gradient operator
    //     Matrix const B = fem::sym_gradient<3>((rhea * Jacobian.inverse()).transpose());
    //
    //     return B.transpose() * D * B * Jacobian.determinant();
    // });
}

Vector femSubmesh::internal_nodal_force(Matrix const& x, int64 const element) const
{
    // auto const& cauchy_stresses = variables(InternalVariables::Tensor::Cauchy);

    auto const[m, n] = std::make_tuple(nodes_per_element(), dofs_per_node());

    Vector fint = Vector::Zero(m * n);

    // sf->quadrature().integrate(Eigen::Map<RowMatrix>(fint.data(), m, n),
    //                            [&](auto const& femval, auto const& l) -> RowMatrix {
    //
    //                                auto const & [N, dN] = femval;
    //
    //                                auto const Jacobian = local_deformation_gradient(dN, x);
    //
    //                                auto const& cauchy_stress = cauchy_stresses[offset(element, l)];
    //
    //                                // Compute the symmetric gradient operator
    //                                auto const Bt = dN * Jacobian.inverse();
    //
    //                                return Bt * cauchy_stress * Jacobian.determinant();
    //                            });
    return fint;
}

std::tuple<List const&, Matrix> femSubmesh::consistent_mass(int64 const element) const
{
    auto const X = material_coordinates->initial_configuration(local_node_list(element));
    return {local_dof_list(element), X};
    // auto const density_0 = cm->intrinsic_material().initial_density();
    //
    // auto m = sf->quadrature().integrate(Matrix::Zero(nodes_per_element(), nodes_per_element()).eval(),
    //                                     [&](auto const& femval, auto const& l) -> Matrix {
    //                                         auto const & [N, dN] = femval;
    //
    //                                         auto const Jacobian = local_deformation_gradient(dN, X);
    //
    //                                         return N * density_0 * N.transpose()
    //                                                * Jacobian.determinant();
    //                                     });
    // return {local_dof_list(element), identity_expansion(m, dofs_per_node())};
}

std::tuple<List const&, Vector> femSubmesh::diagonal_mass(int64 const element) const
{
    auto const & [dofs, consistent_m] = this->consistent_mass(element);

    Vector diagonal_m(consistent_m.rows());
    for (auto i = 0; i < consistent_m.rows(); ++i)
    {
        diagonal_m(i) = consistent_m.row(i).sum();
    }
    return {local_dof_list(element), diagonal_m};
}

void femSubmesh::update_internal_variables(double const time_step_size)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    update_deformation_measures();

    update_Jacobian_determinants();

    cm->update_internal_variables(time_step_size);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }
}

void femSubmesh::update_deformation_measures()
{
    // auto& H_list = variables(InternalVariables::Tensor::DisplacementGradient);
    // auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    //
    // // #pragma omp parallel for
    // for (auto element = 0; element < elements(); ++element)
    // {
    //     // Gather the material coordinates
    //     auto const X = project_to_plane(material_coordinates->initial_configuration(local_node_list(element)));
    //     auto const x = project_to_plane(material_coordinates->current_configuration(local_node_list(element)));
    //
    //     sf->quadrature().for_each([&](auto const& femval, const auto& l) {
    //
    //         auto const & [N, rhea] = femval;
    //
    //         // Local deformation gradient for the initial configuration
    //         Matrix2 const F_0 = local_deformation_gradient(rhea, X);
    //         Matrix2 const F = local_deformation_gradient(rhea, x);
    //
    //         // Gradient operator in index notation
    //         auto const& B_0t = rhea * F_0.inverse();
    //
    //         // Displacement gradient
    //         Matrix2 const H = (x - X) * B_0t;
    //
    //         H_list[offset(element, l)] = H;
    //         F_list[offset(element, l)] = F * F_0.inverse();
    //     });
    // }
}

void femSubmesh::update_Jacobian_determinants()
{
    // auto const& deformation_gradients = variables(InternalVariables::Tensor::DeformationGradient);
    // auto& deformation_gradient_determinants = variables(InternalVariables::Scalar::DetF);
    //
    // deformation_gradient_determinants = deformation_gradients
    //                                     | ranges::view::transform(
    //                                           [](auto const& F) { return F.determinant(); });
    //
    // auto const found = ranges::find_if(deformation_gradient_determinants,
    //                                    [](auto const& detF) { return detF <= 0.0; });
    //
    // if (found != deformation_gradient_determinants.end())
    // {
    //     auto const i = std::distance(deformation_gradient_determinants.begin(), found);
    //
    //     auto const[element, quadrature_point] = std::div(i, sf->quadrature().points());
    //
    //     throw computational_error("Positive Jacobian assumption violated at element "
    //                               + std::to_string(element) + " and local quadrature point "
    //                               + std::to_string(quadrature_point) + " (" + std::to_string(*found)
    //                               + ")");
    // }
}

int64 femSubmesh::offset(int64 const element, int64 const quadrature_point) const
{
    return sf->quadrature().points() * element + quadrature_point;
}
}
