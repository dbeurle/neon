
#include <catch.hpp>

#include "mesh/basic_mesh.hpp"
#include "mesh/mechanical/beam/fem_submesh.hpp"

#include "io/json.hpp"

#include <iostream>

constexpr auto ZERO_MARGIN = 1.0e-5;

std::string simple_beam_element()
{
    // clang-format off
    return "{ \
       \"Elements\" : [ \
          { \
             \"Name\" : \"beam\", \
             \"NodalConnectivity\" : [ \
                [ 0, 2 ], \
                [ 2, 1 ] \
             ], \
             \"Type\" : 1 \
         }, \
          { \
             \"Name\" : \"fix\", \
             \"NodalConnectivity\" : [ \
                [ 0 ] \
             ], \
             \"Type\" : 15 \
         }, \
          { \
             \"Name\" : \"load\", \
             \"NodalConnectivity\" : [ \
                [ 1 ] \
             ], \
             \"Type\" : 15 \
         } \
       ], \
       \"Nodes\" : [ \
          { \
             \"Coordinates\" : [ \
                [ 0, 0, 0 ], \
                [ 1, 0, 0 ], \
                [ 0.49999999999867178, 0, 0 ] \
             ] \
         } \
       ] \
    }";
    // clang-format on
}

neon::matrix compute_linear_reduced_quadrature_stiffness(double const element_length,
                                                         double const area,
                                                         double const shear_area_1,
                                                         double const shear_area_2,
                                                         double const inertia_1,
                                                         double const inertia_2,
                                                         double const shear_modulus,
                                                         double const elastic_modulus)
{
    auto h{element_length};

    auto const J = inertia_1 + inertia_2;

    auto const k1 = shear_modulus * shear_area_1 / h;
    auto const k2 = shear_modulus * shear_area_2 / h;
    auto const k3 = elastic_modulus * area / h;
    auto const k4 = shear_modulus * shear_area_2 * h / 4.0 + elastic_modulus * inertia_1 / h;
    auto const k5 = shear_modulus * shear_area_2 * h / 4.0 - elastic_modulus * inertia_1 / h;
    auto const k6 = shear_modulus * shear_area_1 * h / 4.0 + elastic_modulus * inertia_2 / h;
    auto const k7 = shear_modulus * shear_area_1 * h / 4.0 - elastic_modulus * inertia_2 / h;
    auto const k8 = shear_modulus * J / h;

    neon::matrix ke = neon::matrix::Zero(12, 12);

    ke(0, 0) = k1;
    ke(1, 1) = k2;
    ke(2, 2) = k3;
    ke(3, 3) = k4;
    ke(4, 4) = k6;
    ke(5, 5) = k8;
    ke(6, 6) = k1;
    ke(7, 7) = k2;
    ke(8, 8) = k3;
    ke(9, 9) = k4;
    ke(10, 10) = k6;
    ke(11, 11) = k8;

    ke(0, 4) = ke(4, 0) = h / 2.0 * k1;
    ke(0, 6) = ke(6, 0) = -k1;
    ke(0, 10) = ke(10, 0) = h / 2.0 * k1;

    ke(1, 3) = ke(3, 1) = -h / 2.0 * k2;
    ke(1, 7) = ke(7, 1) = -k2;
    ke(1, 9) = ke(9, 1) = -h / 2.0 * k2;

    ke(2, 8) = ke(8, 2) = -k3;

    ke(3, 7) = ke(7, 3) = h / 2.0 * k2;
    ke(3, 9) = ke(9, 3) = k5;

    ke(4, 6) = ke(6, 4) = -h / 2.0 * k1;
    ke(4, 10) = ke(10, 4) = k7;

    ke(5, 11) = ke(11, 5) = -k8;

    ke(6, 10) = ke(10, 6) = -h / 2.0 * k1;

    ke(7, 9) = ke(9, 7) = h / 2.0 * k2;

    return ke;
}

TEST_CASE("Reduced integration linear element")
{
    // The stiffness matrix can be computed analytically for linear shape functions.
    // Here we will compare the numerically evaluated stiffness matrix against
    // the computational model

    // A test geometry of a square bar is used with base SI units.

    //      a = 1m
    //  <---------->
    //  |----------|  ^
    //  |          |  |   b = 1m
    //  |          |  |
    //  |----------|  \/
    //
    // SECTION PROPERTIES:
    //  I_1 = I_2 = 1/12
    //
    // MATERIAL PROPERTIES:
    // Elastic modulus is 200.0e9 MPa and Poisson's ratio is 0.3 which gives a
    // shear modulus of E / (2 * (1 + nu)) = 7.69e10
    neon::matrix const precomputed_ke = compute_linear_reduced_quadrature_stiffness(0.5,
                                                                                    1.0,
                                                                                    1.0,
                                                                                    1.0,
                                                                                    1.0 / 12.0,
                                                                                    1.0 / 12.0,
                                                                                    200.0e9 / 2.0
                                                                                        / (1.0 + 0.3),
                                                                                    200.0e9);

    std::cout << precomputed_ke << '\n' << '\n';

    neon::basic_mesh beam_mesh(json::parse(simple_beam_element()));
    neon::nodal_coordinates beam_nodes(json::parse(simple_beam_element()));

    auto coordinates = std::make_shared<neon::material_coordinates>(beam_nodes.coordinates());

    REQUIRE(coordinates->size() == 3);

    for (auto const& submesh : beam_mesh.meshes("beam"))
    {
        neon::mechanical::beam::fem_submesh mesh(json::parse("{\"Name\" : \"steel\", \
                                                              \"ElasticModulus\" : 200.0e9, \
                                                              \"PoissonsRatio\" : 0.3 \
                                                              }"),
                                                 json::parse("{\"ElementOptions\" : { \
                                                                 \"Quadrature\" : \"Full\"} \
                                                             }"),
                                                 coordinates,
                                                 submesh);

        mesh.update_internal_variables();

        auto const [dofs, k_e] = mesh.tangent_stiffness(0);

        std::cout << "numerical\n" << k_e << '\n';

        REQUIRE((precomputed_ke - k_e).norm() / (k_e).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
}
