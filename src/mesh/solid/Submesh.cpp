
#include "Submesh.hpp"

#include "PreprocessorExceptions.hpp"
#include "constitutive/Hyperelastic.hpp"
#include "interpolations/Hexahedron8.hpp"
#include "mesh/solid/MaterialCoordinates.hpp"
#include "numeric/Operators.hpp"

#include <chrono>
#include <json/json.h>

#include <range/v3/view.hpp>

namespace neon::solid
{
femSubmesh::femSubmesh(Json::Value const& material_data,
                       Json::Value const& simulation_data,
                       std::shared_ptr<MaterialCoordinates>& material_coordinates,
                       SubMesh const& submesh)
    : neon::SubMesh(submesh),
      material_coordinates(material_coordinates),
      sf(make_shape_function(simulation_data)),
      variables(elements() * sf->quadrature().points()),
      cm(make_constitutive_model(material_data, simulation_data))
{
    // Allocate storage for the displacement gradient
    variables.add(InternalVariables::Tensor::DisplacementGradient,
                  InternalVariables::Tensor::DeformationGradient);

    variables.add(InternalVariables::Scalar::DetJ);

    allocate_dof_list(this->dofs_per_node());
}

std::tuple<List const&, Matrix> femSubmesh::tangent_stiffness(int element) const
{
    auto x = material_coordinates->current_configuration(local_node_list(element));

    Matrix kgeo = geometric_tangent_stiffness(x, element);
    Matrix kmat = material_tangent_stiffness(x, element);

    return {local_dof_list(element), kgeo + kmat};
}

std::tuple<List const&, Vector> femSubmesh::internal_force(int element) const
{
    auto X = material_coordinates->initial_configuration(local_node_list(element));
    auto x = material_coordinates->current_configuration(local_node_list(element));

    return {local_dof_list(element), internal_nodal_force(x, element)};
}

Matrix femSubmesh::geometric_tangent_stiffness(Matrix const& configuration, int element) const
{
    auto const& kirchhoffs = variables[InternalVariables::Tensor::Kirchhoff];

    auto n = nodes_per_element();

    Matrix const kgeo =
        sf->quadrature().integrate(Matrix::Zero(n, n), [&](auto const& femval, auto const& l) -> Matrix {
            auto const & [ N, rhea ] = femval;

            auto const Jacobian = local_deformation_gradient(rhea, configuration);

            auto const j = Jacobian.determinant();

            Matrix3 cauchy_stress = kirchhoffs[offset(element, l)] / j;

            // Compute the symmetric gradient operator
            const auto L = (rhea * Jacobian.inverse()).transpose();

            return L.transpose() * cauchy_stress * L * j;
        });
    return identity_expansion(kgeo, dofs_per_node(), nodes_per_element());
}

Matrix femSubmesh::material_tangent_stiffness(Matrix const& configuration, int element) const
{
    auto const local_dofs = nodes_per_element() * dofs_per_node();

    auto const& D_Vec = variables[InternalVariables::Matrix::MaterialTangent];

    Matrix kmat = Matrix::Zero(local_dofs, local_dofs);

    return sf->quadrature().integrate(kmat, [&](auto const& femval, auto const& l) -> Matrix {

        auto const & [ N, rhea ] = femval;

        auto const& D = D_Vec[offset(element, l)];

        Matrix3 const Jacobian = local_deformation_gradient(rhea, configuration);

        auto const detJ = Jacobian.determinant();

        // Compute the symmetric gradient operator
        Matrix const B = fem::SymGradient<3>((rhea * Jacobian.inverse()).transpose());

        return B.transpose() * D * B * detJ;
    });
}

Vector femSubmesh::internal_nodal_force(Matrix const& configuration, int element) const
{
    using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    auto const& kirchhoff = variables[InternalVariables::Tensor::Kirchhoff];

    auto const n = nodes_per_element();
    auto const m = dofs_per_node();

    RowMatrix fint =
        sf->quadrature().integrate(RowMatrix::Zero(n, m),
                                   [&](auto const& femval, auto const& l) -> RowMatrix {
                                       auto const & [ N, dN ] = femval;

                                       auto const Jacobian =
                                           local_deformation_gradient(dN, configuration);

                                       auto const j = Jacobian.determinant();

                                       auto const σ = kirchhoff[offset(element, l)] / j;

                                       // Compute the symmetric gradient operator
                                       const auto Bt = dN * Jacobian.inverse();

                                       return Bt * σ * j;
                                   });
    return Eigen::Map<RowMatrix>(fint.data(), n * m, 1);
}

void femSubmesh::update_internal_variables()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Updating the internal variables..." << std::flush;

    update_deformation_measures();
    update_Jacobian_determinants();
    check_element_distortion();

    cm->update_internal_variables();
    cm->update_continuum_tangent();

    variables.commit();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "took " << elapsed_seconds.count() << "s\n";
}

void femSubmesh::update_deformation_measures()
{
    // FIXME These is a hack because of the decomposition and lambda functions
    auto& Hvec = variables(InternalVariables::Tensor::DisplacementGradient);
    auto& Fvec = variables(InternalVariables::Tensor::DeformationGradient);

    // #pragma omp parallel for
    for (auto element = 0; element < elements(); ++element)
    {
        // Gather the material coordinates
        auto const X = material_coordinates->initial_configuration(local_node_list(element));
        auto const x = material_coordinates->current_configuration(local_node_list(element));

        sf->quadrature().for_each([&](auto const& femval, const auto& l) {

            auto const & [ N, rhea ] = femval;

            // Local deformation gradient for the initial configuration
            Matrix3 const& F_0 = local_deformation_gradient(rhea, X);

            // Gradient operator in index notation
            auto const& B_0t = rhea * F_0.inverse();

            // Nodal displacement matrix
            auto const u = x - X;

            // Strain operator
            Matrix3 const H = u * B_0t;

            Hvec[offset(element, l)] = H;
            Fvec[offset(element, l)] = Matrix3::Identity() + H;
        });
    }
}

void femSubmesh::update_deformation_gradient(int element, Matrix const& X, Matrix const& x)
{
    auto& Fvec = variables(InternalVariables::Tensor::DeformationGradient);

    sf->quadrature().for_each([&](auto const& femval, const auto& l) {

        auto const & [ N, rhea ] = femval;

        Fvec[offset(element, l)] = deformation_gradient(rhea, x, X);
    });
}

void femSubmesh::update_Jacobian_determinants()
{
    auto const& Fvec = variables(InternalVariables::Tensor::DeformationGradient);

    auto& detJvec = variables(InternalVariables::Scalar::DetJ);

    detJvec = Fvec | ranges::view::transform([](auto const& F) { return F.determinant(); });
}

void femSubmesh::check_element_distortion() const
{
    auto const& detJvec = variables(InternalVariables::Scalar::DetJ);

    auto found = ranges::find_if(detJvec, [](auto const& detJ) { return detJ < 0.0; });

    if (found != detJvec.end())
    {
        auto const i = std::distance(detJvec.begin(), found);

        auto const element = std::floor(static_cast<double>(i) / sf->quadrature().points());
        auto const quadrature_point = i % sf->quadrature().points();
        throw DistortedElement(element, quadrature_point);
    }
}

void femSubmesh::allocate_dof_list(int const nodal_dofs)
{
    using namespace ranges;

    dof_list = nodal_connectivity | view::transform([=](auto const& node_list) {
                   return view::for_each(node_list, [=](int const local_node) {
                       return view::ints(0, nodal_dofs) | view::transform([=](int const nodal_dof) {
                                  return local_node * nodal_dofs + nodal_dof;
                              });
                   });
               });
}

std::unique_ptr<ConstitutiveModel> femSubmesh::make_constitutive_model(
    Json::Value const& material_data,
    Json::Value const& simulation_data)
{
    if (simulation_data["ConstitutiveModel"].empty())
    {
        std::cout << simulation_data << std::endl;
        throw EmptyFieldException("Part: ConstitutiveModel");
    }

    auto const& model_name = simulation_data["ConstitutiveModel"].asString();

    if (model_name == "NeoHooke")
    {
        return std::make_unique<NeoHooke>(variables, material_data);
    }

    std::cerr << "The model name " << model_name << " is not recognised\n";
    std::abort();
    return nullptr;
}

std::unique_ptr<VolumeInterpolation> femSubmesh::make_shape_function(Json::Value const& simulation_data)
{
    std::cout << "Allocating a volume interpolation function" << std::endl;

    if (!simulation_data.isMember("ElementOptions"))
    {
        throw EmptyFieldException("Part: ElementOptions");
    }

    auto is_reduced = simulation_data["ElementOptions"]["Quadrature"].empty()
                          ? false
                          : simulation_data["ElementOptions"]["Quadrature"].asString() == "Reduced";

    switch (topology())
    {
        case ElementTopology::Hexahedron8:
        {
            return std::make_unique<Hexahedron8>(is_reduced ? HexahedronQuadrature::Rule::OnePoint
                                                            : HexahedronQuadrature::Rule::EightPoint);
            break;
        }
        break;
        case ElementTopology::Tetrahedron4:
            break;
        case ElementTopology::Prism6:
        case ElementTopology::Pyramid5:
        case ElementTopology::Tetrahedron20:
        case ElementTopology::Tetrahedron35:
        case ElementTopology::Tetrahedron56:
        case ElementTopology::Hexahedron64:
        case ElementTopology::Hexahedron125:
        case ElementTopology::Tetrahedron10:
        case ElementTopology::Hexahedron27:
        case ElementTopology::Prism18:
        case ElementTopology::Pyramid14:
        case ElementTopology::Hexahedron20:
        case ElementTopology::Prism15:
        case ElementTopology::Pyramid13:
            std::cerr << "Element shape not implemented inside the solid femSubmesh\n";
            std::terminate();
            break;
        case ElementTopology::Invalid:
        case ElementTopology::Point:
        case ElementTopology::Line2:
        case ElementTopology::Line3:
        case ElementTopology::Triangle3:
        case ElementTopology::Triangle6:
        case ElementTopology::Triangle9:
        case ElementTopology::Quadrilateral4:
        case ElementTopology::Quadrilateral8:
        case ElementTopology::Quadrilateral9:
        case ElementTopology::Triangle10:
        case ElementTopology::Triangle12:
        case ElementTopology::Triangle15:
        case ElementTopology::Triangle15_IC:
        case ElementTopology::Triangle21:
        case ElementTopology::Edge4:
        case ElementTopology::Edge5:
        case ElementTopology::Edge6:
        default:
            break;
    }
    return nullptr;
}

int femSubmesh::offset(int element, int quadrature_point) const
{
    return sf->quadrature().points() * element + quadrature_point;
}
}
