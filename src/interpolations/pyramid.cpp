
#include "pyramid.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
pyramid5::pyramid5(pyramid_quadrature::point const p)
    : volume_interpolation(std::make_unique<pyramid_quadrature>(p), 5)
{
    this->precompute_shape_functions();
}

void pyramid5::precompute_shape_functions()
{
    using coordinate_type = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinate_type, 5> constexpr local_coordinates{{{0, 1.0, 0.0, -1.0},
                                                                {1, 0.0, 1.0, -1.0},
                                                                {2, 0.0, 0.0, -1.0},
                                                                {3, 1.0, 0.0, 1.0},
                                                                {4, 0.0, 1.0, 1.0}}};

    matrix N_matrix(m_quadrature->points(), number_of_nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 4);

    m_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, u1, u2, u3] = coordinate;

        vector N(5);
        matrix dN(5, 3);

        N(0) = -(u1 - u3 + 1) * (u2 - u3 + 1) / (4 * (u3 - 1));
        N(1) = (u1 + u3 - 1) * (u2 - u3 + 1) / (4 * (u3 - 1));
        N(2) = -(u1 + u3 - 1) * (u2 + u3 - 1) / (4 * (u3 - 1));
        N(3) = (u1 - u3 + 1) * (u2 + u3 - 1) / (4 * (u3 - 1));
        N(4) = u3;

        dN(0, 0) = -(u2 - u3 + 1) / (4 * u3 - 4);
        dN(1, 0) = (u2 - u3 + 1) / (4 * (u3 - 1));
        dN(2, 0) = -(u2 + u3 - 1) / (4 * u3 - 4);
        dN(3, 0) = (u2 + u3 - 1) / (4 * (u3 - 1));
        dN(4, 0) = 0;

        dN(0, 1) = -(u1 - u3 + 1) / (4 * u3 - 4);
        dN(1, 1) = (u1 + u3 - 1) / (4 * (u3 - 1));
        dN(2, 1) = -(u1 + u3 - 1) / (4 * u3 - 4);
        dN(3, 1) = (u1 - u3 + 1) / (4 * (u3 - 1));
        dN(4, 1) = 0;

        dN(0, 2) = (u1 * u2 - std::pow(u3, 2) + 2 * u3 - 1) / (4 * (std::pow(u3, 2) - 2 * u3 + 1));
        dN(1, 2) = -(u1 * u2 + std::pow(u3, 2) - 2 * u3 + 1) / (4 * std::pow(u3, 2) - 8 * u3 + 4);
        dN(2, 2) = (u1 * u2 - std::pow(u3, 2) + 2 * u3 - 1) / (4 * (std::pow(u3, 2) - 2 * u3 + 1));
        dN(3, 2) = -(u1 * u2 + std::pow(u3, 2) - 2 * u3 + 1) / (4 * std::pow(u3, 2) - 8 * u3 + 4);
        dN(4, 2) = 1;

        local_quadrature_coordinates(l, 0) = u1;
        local_quadrature_coordinates(l, 1) = u2;
        local_quadrature_coordinates(l, 2) = u3;

        N_matrix.row(l) = N;

        return std::make_tuple(N, dN);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(number_of_nodes(), 4);

    for (auto const& [a, r_a, s_a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r_a;
        local_nodal_coordinates(a, 1) = s_a;
        local_nodal_coordinates(a, 2) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double pyramid5::compute_measure(matrix3x const& nodal_coordinates) const
{
    return m_quadrature->integrate(0.0, [&](auto const& femval, auto) {
        auto const& [N, dN] = femval;

        matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}

pyramid13::pyramid13(pyramid_quadrature::point const p)
    : volume_interpolation(std::make_unique<pyramid_quadrature>(p), 13)
{
    this->precompute_shape_functions();
}

void pyramid13::precompute_shape_functions()
{
    using coordinate_type = std::tuple<int, double, double, double>;

    // Initialise nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinate_type, 13> constexpr local_coordinates{{{0, 1.0, 0.0, -1.0},
                                                                 {1, 0.0, 1.0, -1.0},
                                                                 {2, 0.0, 0.0, -1.0},
                                                                 {3, 0.5, 0.5, -1.0},
                                                                 {4, 0.0, 0.5, -1.0},
                                                                 {5, 0.5, 0.0, -1.0},
                                                                 //
                                                                 {6, 1.0, 0.0, 0.0},
                                                                 {7, 0.0, 1.0, 0.0},
                                                                 {8, 0.0, 0.0, 0.0},
                                                                 //
                                                                 {9, 1.0, 0.0, 1.0},
                                                                 {10, 0.0, 1.0, 1.0},
                                                                 {11, 0.0, 0.0, 1.0},
                                                                 {12, 0.0, 0.0, 1.0}}};

    matrix N_matrix(m_quadrature->points(), number_of_nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 4);

    m_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, u1, u2, u3] = coordinate;

        vector N(13);
        matrix dN(13, 3);

        N(0) = -(u1 + u2 - 1) * (u1 - u3 + 1) * (u2 - u3 + 1) / (4 * (u3 - 1));
        N(1) = (u1 - u3 + 1) * (u1 + u3 - 1) * (u2 - u3 + 1) / (2 * (u3 - 1));
        N(2) = -(u1 - u2 + 1) * (u1 + u3 - 1) * (u2 - u3 + 1) / (4 * (u3 - 1));
        N(3) = -(u1 + u3 - 1) * (u2 - u3 + 1) * (u2 + u3 - 1) / (2 * (u3 - 1));
        N(4) = (u1 + u2 + 1) * (u1 + u3 - 1) * (u2 + u3 - 1) / (4 * (u3 - 1));
        N(5) = -(u1 - u3 + 1) * (u1 + u3 - 1) * (u2 + u3 - 1) / (2 * (u3 - 1));
        N(6) = (u1 - u2 - 1) * (u1 - u3 + 1) * (u2 + u3 - 1) / (4 * (u3 - 1));
        N(7) = (u1 - u3 + 1) * (u2 - u3 + 1) * (u2 + u3 - 1) / (2 * (u3 - 1));
        N(8) = -u3 * (u1 - u3 + 1) * (u2 - u3 + 1) / (u3 - 1);
        N(9) = u3 * (u1 + u3 - 1) * (u2 - u3 + 1) / (u3 - 1);
        N(10) = -u3 * (u1 + u3 - 1) * (u2 + u3 - 1) / (u3 - 1);
        N(11) = u3 * (u1 - u3 + 1) * (u2 + u3 - 1) / (u3 - 1);
        N(12) = u3 * (2 * u3 - 1);

        dN(0, 0) = -(2 * u1 * u2 - 2 * u1 * u3 + 2 * u1 + std::pow(u2, 2) - 2 * u2 * u3 + u2
                     + std::pow(u3, 2) - u3)
                   / (4 * u3 - 4);
        dN(1, 0) = u1 * (u2 - u3 + 1) / (u3 - 1);
        dN(2, 0) = -(2 * u1 * u2 - 2 * u1 * u3 + 2 * u1 - std::pow(u2, 2) + 2 * u2 * u3 - u2
                     - std::pow(u3, 2) + u3)
                   / (4 * u3 - 4);
        dN(3, 0) = -(u2 - u3 + 1) * (u2 + u3 - 1) / (2 * u3 - 2);
        dN(4, 0) = (2 * u1 * u2 + 2 * u1 * u3 - 2 * u1 + std::pow(u2, 2) + 2 * u2 * u3 - u2
                    + std::pow(u3, 2) - u3)
                   / (4 * (u3 - 1));
        dN(5, 0) = -u1 * (u2 + u3 - 1) / (u3 - 1);
        dN(6, 0) = (2 * u1 * u2 + 2 * u1 * u3 - 2 * u1 - std::pow(u2, 2) - 2 * u2 * u3 + u2
                    - std::pow(u3, 2) + u3)
                   / (4 * (u3 - 1));
        dN(7, 0) = (u2 - u3 + 1) * (u2 + u3 - 1) / (2 * (u3 - 1));
        dN(8, 0) = -u3 * (u2 - u3 + 1) / (u3 - 1);
        dN(9, 0) = u3 * (u2 - u3 + 1) / (u3 - 1);
        dN(10, 0) = -u3 * (u2 + u3 - 1) / (u3 - 1);
        dN(11, 0) = u3 * (u2 + u3 - 1) / (u3 - 1);
        dN(12, 0) = 0.0;

        dN(0, 1) = -(std::pow(u1, 2) + 2 * u1 * u2 - 2 * u1 * u3 + u1 - 2 * u2 * u3 + 2 * u2
                     + std::pow(u3, 2) - u3)
                   / (4 * u3 - 4);
        dN(1, 1) = (u1 - u3 + 1) * (u1 + u3 - 1) / (2 * (u3 - 1));
        dN(2, 1) = -(std::pow(u1, 2) - 2 * u1 * u2 + 2 * u1 * u3 - u1 - 2 * u2 * u3 + 2 * u2
                     + std::pow(u3, 2) - u3)
                   / (4 * u3 - 4);
        dN(3, 1) = -u2 * (u1 + u3 - 1) / (u3 - 1);
        dN(4, 1) = (std::pow(u1, 2) + 2 * u1 * u2 + 2 * u1 * u3 - u1 + 2 * u2 * u3 - 2 * u2
                    + std::pow(u3, 2) - u3)
                   / (4 * (u3 - 1));
        dN(5, 1) = -(u1 - u3 + 1) * (u1 + u3 - 1) / (2 * u3 - 2);
        dN(6, 1) = (std::pow(u1, 2) - 2 * u1 * u2 - 2 * u1 * u3 + u1 + 2 * u2 * u3 - 2 * u2
                    + std::pow(u3, 2) - u3)
                   / (4 * (u3 - 1));
        dN(7, 1) = u2 * (u1 - u3 + 1) / (u3 - 1);
        dN(8, 1) = -u3 * (u1 - u3 + 1) / (u3 - 1);
        dN(9, 1) = u3 * (u1 + u3 - 1) / (u3 - 1);
        dN(10, 1) = -u3 * (u1 + u3 - 1) / (u3 - 1);
        dN(11, 1) = u3 * (u1 - u3 + 1) / (u3 - 1);
        dN(12, 1) = 0.0;

        dN(0, 2) = -(u1 + u2 - 1) * (-u1 * u2 * u3 + u1 * u2 * (u3 - 1) + std::pow(u3 - 1, 2))
                   / (4 * std::pow(u3 - 1, 2));
        dN(1, 2) = (-std::pow(u1, 2) * u2 / 2 - u2 * std::pow(u3, 2) / 2 + u2 * u3 - u2 / 2
                    + std::pow(u3, 3) - 3 * std::pow(u3, 2) + 3 * u3 - 1)
                   / (std::pow(u3, 2) - 2 * u3 + 1);
        dN(2, 2) = (u1 - u2 + 1) * (u1 * u2 * u3 - u1 * u2 * (u3 - 1) + std::pow(u3 - 1, 2))
                   / (4 * std::pow(u3 - 1, 2));
        dN(3, 2) = (u1 * std::pow(u2, 2) + u1 * std::pow(u3, 2) - 2 * u1 * u3 + u1
                    + 2 * std::pow(u3, 3) - 6 * std::pow(u3, 2) + 6 * u3 - 2)
                   / (2 * (std::pow(u3, 2) - 2 * u3 + 1));
        dN(4, 2) = (u1 + u2 + 1) * (-u1 * u2 * u3 + u1 * u2 * (u3 - 1) + std::pow(u3 - 1, 2))
                   / (4 * std::pow(u3 - 1, 2));
        dN(5, 2) = (std::pow(u1, 2) * u2 + u2 * std::pow(u3, 2) - 2 * u2 * u3 + u2
                    + 2 * std::pow(u3, 3) - 6 * std::pow(u3, 2) + 6 * u3 - 2)
                   / (2 * (std::pow(u3, 2) - 2 * u3 + 1));
        dN(6, 2) = (-u1 + u2 + 1) * (u1 * u2 * u3 - u1 * u2 * (u3 - 1) + std::pow(u3 - 1, 2))
                   / (4 * std::pow(u3 - 1, 2));
        dN(7, 2) = (-u1 * std::pow(u2, 2) / 2 - u1 * std::pow(u3, 2) / 2 + u1 * u3 - u1 / 2
                    + std::pow(u3, 3) - 3 * std::pow(u3, 2) + 3 * u3 - 1)
                   / (std::pow(u3, 2) - 2 * u3 + 1);
        dN(8, 2) = (u3 * (u1 - u3 + 1) * (u2 - u3 + 1)
                    + (u3 - 1)
                          * (u3 * (u1 - u3 + 1) + u3 * (u2 - u3 + 1) - (u1 - u3 + 1) * (u2 - u3 + 1)))
                   / std::pow(u3 - 1, 2);
        dN(9,
           2) = (-u3 * (u1 + u3 - 1) * (u2 - u3 + 1)
                 + (u3 - 1)
                       * (-u3 * (u1 + u3 - 1) + u3 * (u2 - u3 + 1) + (u1 + u3 - 1) * (u2 - u3 + 1)))
                / std::pow(u3 - 1, 2);
        dN(10,
           2) = (u3 * (u1 + u3 - 1) * (u2 + u3 - 1)
                 - (u3 - 1)
                       * (u3 * (u1 + u3 - 1) + u3 * (u2 + u3 - 1) + (u1 + u3 - 1) * (u2 + u3 - 1)))
                / std::pow(u3 - 1, 2);
        dN(11,
           2) = (-u3 * (u1 - u3 + 1) * (u2 + u3 - 1)
                 + (u3 - 1)
                       * (u3 * (u1 - u3 + 1) - u3 * (u2 + u3 - 1) + (u1 - u3 + 1) * (u2 + u3 - 1)))
                / std::pow(u3 - 1, 2);

        dN(12, 2) = 4 * u3 - 1;

        local_quadrature_coordinates(l, 0) = u1;
        local_quadrature_coordinates(l, 1) = u2;
        local_quadrature_coordinates(l, 2) = u3;

        N_matrix.row(l) = N;

        return std::make_tuple(N, dN);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(number_of_nodes(), 4);

    for (auto const& [a, r_a, s_a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r_a;
        local_nodal_coordinates(a, 1) = s_a;
        local_nodal_coordinates(a, 2) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double pyramid13::compute_measure(matrix3x const& nodal_coordinates) const
{
    return m_quadrature->integrate(0.0, [&](auto const& femval, auto) {
        auto const& [N, dN] = femval;

        matrix3 const jacobian = nodal_coordinates * dN;

        return jacobian.determinant();
    });
}
}
