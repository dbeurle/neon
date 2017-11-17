
#pragma once

#include "mesh/generic/Neumann.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/diffusion/Submesh.hpp"

namespace neon::diffusion
{
/**
 * SurfaceFlux is a group of the same elements on the same boundary. These
 * elements are responsible for computing their elemental right hand side contributions
 * with the corresponding shape function.  These are required to be stored in a parent
 * container with the other groups from the collective boundary \sa SurfaceFluxBoundary
 */
using SurfaceFlux = SurfaceLoad<SurfaceInterpolation>;

/**
 * NewtonConvection is responsible for computing the additional term to the
 * stiffness matrix and right hand side for modelling physics such as the Newton's
 * law of cooling given by
 *
   \f{align*}{
     h &= \lambda u - q_i n_i
   \f}
 * where \f$ h \f$ is the heat flux (\f$ W/(m^2) \f$), \f$ \lambda \f$ is the heat
 * transfer coefficient \f$ W / (m^2 K) \f$ and \f$ q_i n_i \f$ is the heat
 * leaving the surface.  The surrounding temperature \f$T_\infty\f$ is computed by
 * \f$ T_\infty = h / \lambda \f$
 */
class NewtonConvection : public SurfaceFlux
{
public:
    explicit NewtonConvection(std::unique_ptr<SurfaceInterpolation>&& sf,
                              std::vector<List> const& nodal_connectivity,
                              std::vector<List> const& dof_list,
                              std::shared_ptr<MaterialCoordinates>& material_coordinates,
                              Json::Value const& time_history,
                              Json::Value const& heat_transfer_coefficients,
                              Json::Value const& prescribed_heat_fluxes);

    /**
     * Compute the element stiffness matrix contributing to the mixed boundary condition
       \f{align*}{
         k_{ab} &= \int_{\Gamma} N_a \lambda N_b d\Gamma
       \f}
     */
    [[nodiscard]] std::tuple<List const&, Matrix> external_stiffness(int const element,
                                                                     double const load_factor) const
    {
        // clang-format off
        auto const X = sf->project_to_plane(material_coordinates->initial_configuration(nodal_connectivity[element]));

        auto const coefficient = interpolate_prescribed_load(coefficient_time_data, load_factor);

        // Perform the computation of the external element stiffness matrix
        auto const k_ext = sf->quadrature().integrate(Matrix::Zero(X.cols(), X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> Matrix {
                                                          auto const & [ N, dN ] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return coefficient * N * N.transpose() * j;
                                                      });
        return {dof_list[element], k_ext};
        // clang-format on
    }

protected:
    std::vector<std::pair<double, double>> coefficient_time_data;
};

/* SurfaceFluxBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 */
class SurfaceFluxBoundary
{
public:
    explicit SurfaceFluxBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                 std::vector<SubMesh> const& submeshes,
                                 Json::Value const& times,
                                 Json::Value const& loads,
                                 Json::Value const& simulation_data)
    {
        for (auto const& mesh : submeshes)
        {
            surface_loads.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                                       mesh.connectivities(),
                                       mesh.connectivities(),
                                       material_coordinates,
                                       times,
                                       loads);
        }
    }

    auto const& boundaries() const { return surface_loads; }

protected:
    std::vector<SurfaceFlux> surface_loads;
};
}
