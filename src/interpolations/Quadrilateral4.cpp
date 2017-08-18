
#include "Quadrilateral4.hpp"

#include <array>
#include <tuple>

namespace neon
{
Quadrilateral4::Quadrilateral4(QuadrilateralQuadrature::Rule rule)
    : SurfaceInterpolation(std::make_unique<QuadrilateralQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Quadrilateral4::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    std::array<NodalCoordinate, 4> constexpr local_coordinates = {
        {{0, -1.0, -1.0}, {1, 1.0, -1.0}, {2, 1.0, 1.0}, {3, -1.0, 1.0}}};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinate) {

        auto const & [ l, xi, eta ] = coordinate;

        Vector N(4);
        Matrix rhea(4, 2);

        for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
        {
            N(a) = 0.25 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta);
            rhea(a, 0) = 0.25 * (1.0 + eta_a * eta) * xi_a;
            rhea(a, 1) = 0.25 * (1.0 + xi_a * xi) * eta_a;
        }

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    Matrix local_nodal_coordinates = Matrix::Ones(nodes(), 3);

    for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double Quadrilateral4::compute_measure(Matrix const& nodal_coordinates)
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const & [ N, dN ] = femval;

        Matrix2 const Jacobian = project_to_plane(nodal_coordinates) * dN;

        return Jacobian.determinant();
    });
}
}
