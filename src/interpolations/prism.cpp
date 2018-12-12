
#include "prism.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
prism6::prism6(prism_quadrature::point const p)
    : volume_interpolation(std::make_unique<prism_quadrature>(p))
{
    this->precompute_shape_functions();
}

void prism6::precompute_shape_functions()
{
    using coordinate_type = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinate_type, 6> constexpr local_coordinates{{{0, 1.0, 0.0, -1.0},
                                                                {1, 0.0, 1.0, -1.0},
                                                                {2, 0.0, 0.0, -1.0},
                                                                {3, 1.0, 0.0, 1.0},
                                                                {4, 0.0, 1.0, 1.0},
                                                                {5, 0.0, 0.0, 1.0}}};

    matrix N_matrix(m_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 4);

    m_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s, xi] = coordinate;

        vector N(6);
        matrix rhea(6, 3);

        auto const t = 1.0 - r - s;

        N(1) = 0.5 * (1.0 - xi) * r;
        N(0) = 0.5 * (1.0 - xi) * s;
        N(2) = 0.5 * (1.0 - xi) * t;
        N(3) = 0.5 * (1.0 + xi) * r;
        N(4) = 0.5 * (1.0 + xi) * s;
        N(5) = 0.5 * (1.0 + xi) * t;

        rhea(0, 0) = 0.5 * (1.0 - xi);
        rhea(1, 0) = 0.0;
        rhea(2, 0) = -0.5 * (1.0 - xi);
        rhea(3, 0) = 0.5 * (1.0 + xi);
        rhea(4, 0) = 0.0;
        rhea(5, 0) = -0.5 * (1.0 + xi);

        rhea(0, 1) = 0.0;
        rhea(1, 1) = 0.5 * (1.0 - xi);
        rhea(2, 1) = -0.5 * (1.0 - xi);
        rhea(3, 1) = 0.0;
        rhea(4, 1) = 0.5 * (1.0 + xi);
        rhea(5, 1) = -0.5 * (1.0 + xi);

        rhea(0, 2) = -0.5 * r;
        rhea(1, 2) = -0.5 * s;
        rhea(2, 2) = -0.5 * t;
        rhea(3, 2) = 0.5 * r;
        rhea(4, 2) = 0.5 * s;
        rhea(5, 2) = 0.5 * t;

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;
        local_quadrature_coordinates(l, 2) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 4);

    for (auto const& [a, r_a, s_a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r_a;
        local_nodal_coordinates(a, 1) = s_a;
        local_nodal_coordinates(a, 2) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double prism6::compute_measure(matrix3x const& nodal_coordinates) const
{
    return m_quadrature->integrate(0.0, [&](auto const& femval, auto) {
        auto const& [N, dN] = femval;

        matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}

prism15::prism15(prism_quadrature::point const p)
    : volume_interpolation(std::make_unique<prism_quadrature>(p))
{
    this->precompute_shape_functions();
}

void prism15::precompute_shape_functions()
{
    using coordinate_type = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinate_type, 15> constexpr local_coordinates{{{0, 1.0, 0.0, -1.0},
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
                                                                 {12, 0.5, 0.5, 1.0},
                                                                 {13, 0.0, 0.5, 1.0},
                                                                 {14, 0.5, 0.0, 1.0}}};

    matrix N_matrix(m_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 4);

    m_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s, xi] = coordinate;

        vector N(15);
        matrix rhea(15, 3);

        N(0) = r * xi * (2 * r - 1) * (xi - 1) / 2.0;
        N(1) = s * xi * (2 * s - 1) * (xi - 1) / 2.0;
        N(2) = xi * (xi - 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
        N(3) = 2 * r * s * xi * (xi - 1);
        N(4) = -2 * s * xi * (xi - 1) * (r + s - 1);
        N(5) = -2 * r * xi * (xi - 1) * (r + s - 1);
        N(6) = -r * (xi - 1) * (xi + 1);
        N(7) = -s * (xi - 1) * (xi + 1);
        N(8) = (xi - 1) * (xi + 1) * (r + s - 1);
        N(9) = r * xi * (2 * r - 1) * (xi + 1) / 2.0;
        N(10) = s * xi * (2 * s - 1) * (xi + 1) / 2.0;
        N(11) = xi * (xi + 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
        N(12) = 2 * r * s * xi * (xi + 1);
        N(13) = -2 * s * xi * (xi + 1) * (r + s - 1);
        N(14) = -2 * r * xi * (xi + 1) * (r + s - 1);

        rhea(0, 0) = xi * (4 * r - 1) * (xi - 1) / 2.0;
        rhea(1, 0) = 0;
        rhea(2, 0) = xi * (xi - 1) * (4 * r + 4 * s - 3) / 2.0;
        rhea(3, 0) = 2 * s * xi * (xi - 1);
        rhea(4, 0) = -2 * s * xi * (xi - 1);
        rhea(5, 0) = -2 * xi * (xi - 1) * (2 * r + s - 1);
        rhea(6, 0) = -(xi - 1) * (xi + 1);
        rhea(7, 0) = 0;
        rhea(8, 0) = (xi - 1) * (xi + 1);
        rhea(9, 0) = xi * (4 * r - 1) * (xi + 1) / 2.0;
        rhea(10, 0) = 0;
        rhea(11, 0) = xi * (xi + 1) * (4 * r + 4 * s - 3) / 2.0;
        rhea(12, 0) = 2 * s * xi * (xi + 1);
        rhea(13, 0) = -2 * s * xi * (xi + 1);
        rhea(14, 0) = -2 * xi * (xi + 1) * (2 * r + s - 1);

        rhea(0, 1) = 0;
        rhea(1, 1) = xi * (4 * s - 1) * (xi - 1) / 2.0;
        rhea(2, 1) = xi * (xi - 1) * (4 * r + 4 * s - 3) / 2.0;
        rhea(3, 1) = 2 * r * xi * (xi - 1);
        rhea(4, 1) = -2 * xi * (xi - 1) * (r + 2 * s - 1);
        rhea(5, 1) = -2 * r * xi * (xi - 1);
        rhea(6, 1) = 0;
        rhea(7, 1) = -(xi - 1) * (xi + 1);
        rhea(8, 1) = (xi - 1) * (xi + 1);
        rhea(9, 1) = 0;
        rhea(10, 1) = xi * (4 * s - 1) * (xi + 1) / 2.0;
        rhea(11, 1) = xi * (xi + 1) * (4 * r + 4 * s - 3) / 2.0;
        rhea(12, 1) = 2 * r * xi * (xi + 1);
        rhea(13, 1) = -2 * xi * (xi + 1) * (r + 2 * s - 1);
        rhea(14, 1) = -2 * r * xi * (xi + 1);

        rhea(0, 2) = r * (2 * r - 1) * (2 * xi - 1) / 2.0;
        rhea(1, 2) = s * (2 * s - 1) * (2 * xi - 1) / 2.0;
        rhea(2, 2) = (2 * xi - 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
        rhea(3, 2) = 2 * r * s * (2 * xi - 1);
        rhea(4, 2) = -2 * s * (2 * xi - 1) * (r + s - 1);
        rhea(5, 2) = -2 * r * (2 * xi - 1) * (r + s - 1);
        rhea(6, 2) = -2 * r * xi;
        rhea(7, 2) = -2 * s * xi;
        rhea(8, 2) = 2 * xi * (r + s - 1);
        rhea(9, 2) = r * (2 * r - 1) * (2 * xi + 1) / 2.0;
        rhea(10, 2) = s * (2 * s - 1) * (2 * xi + 1) / 2.0;
        rhea(11, 2) = (2 * xi + 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
        rhea(12, 2) = 2 * r * s * (2 * xi + 1);
        rhea(13, 2) = -2 * s * (2 * xi + 1) * (r + s - 1);
        rhea(14, 2) = -2 * r * (2 * xi + 1) * (r + s - 1);

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;
        local_quadrature_coordinates(l, 2) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 4);

    for (auto const& [a, r_a, s_a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r_a;
        local_nodal_coordinates(a, 1) = s_a;
        local_nodal_coordinates(a, 2) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double prism15::compute_measure(matrix3x const& nodal_coordinates) const
{
    return m_quadrature->integrate(0.0, [&](auto const& femval, auto) {
        auto const& [N, dN] = femval;

        matrix3 const jacobian = nodal_coordinates * dN;

        return jacobian.determinant();
    });
}
}
