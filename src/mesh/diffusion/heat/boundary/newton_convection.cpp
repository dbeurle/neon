
#include "newton_convection.hpp"

#include "math/jacobian_determinant.hpp"
#include "io/json.hpp"

#include <utility>

namespace neon::diffusion::boundary
{
newton_convection::newton_convection(indices node_indices,
                                     indices dof_indices,
                                     std::shared_ptr<material_coordinates>& coordinates,
                                     json const& times,
                                     json const& heat_fluxes,
                                     json const& heat_transfer_coefficients,
                                     element_topology const topology,
                                     json const& element_data)
    : heat_flux(std::move(node_indices),
                std::move(dof_indices),
                coordinates,
                times,
                heat_fluxes,
                topology,
                element_data),
      m_heat_transfer_coefficients{heat_transfer_coefficients.get<std::vector<double>>()}
{
    if (m_heat_transfer_coefficients.size() != m_times.size())
    {
        throw std::domain_error("Heat transfer coefficients need to be the same size");
    }
}

auto newton_convection::external_stiffness(std::int64_t const element,
                                           double const load_factor) const -> matrix const&
{
    auto const configuration = coordinates->initial_configuration(local_node_view(element));

    thread_local matrix k_ext;

    k_ext.resize(configuration.cols(), configuration.cols());

    // Perform the computation of the external element stiffness matrix
    linear_form.integrate(k_ext.setZero(), [&](auto const& values, auto) -> matrix {
        auto const& [N, dN] = values;

        auto const determinant = jacobian_determinant(configuration * dN);

        return N * N.transpose() * determinant;
    });

    k_ext *= interpolate(m_times, m_heat_transfer_coefficients, load_factor);

    return k_ext;
}
}
