
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/basic_mesh.hpp"
#include "mesh/element_topology.hpp"
#include "mesh/node_ordering_adapter.hpp"

#include "mesh/material_coordinates.hpp"
#include "mesh/mechanical/solid/fem_mesh.hpp"
#include "mesh/mechanical/solid/fem_submesh.hpp"

#include "io/json.hpp"

#include <range/v3/view.hpp>

#include <sstream>

#include "CubeJson.hpp"

using namespace neon;
using namespace ranges;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Basic mesh test")
{
    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    basic_mesh basic_mesh(json::parse(cube_mesh_json()));

    nodal_coordinates nodal_coordinates(json::parse(cube_mesh_json()));

    auto constexpr number_of_nodes = 64;

    REQUIRE(nodal_coordinates.size() == number_of_nodes);

    SECTION("Test corner vertices")
    {
        REQUIRE((nodal_coordinates.coordinates(0) - vector3(0.0, 0.0, 0.0)).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((nodal_coordinates.coordinates(1) - vector3(1.0, 0.0, 0.0)).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((nodal_coordinates.coordinates(2) - vector3(0.0, 1.0, 0.0)).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((nodal_coordinates.coordinates(3) - vector3(1.0, 1.0, 0.0)).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Test mesh data for boundary and volume elements")
    {
        for (auto const& mesh : basic_mesh.meshes("bottom"))
        {
            REQUIRE(mesh.topology() == element_topology::quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 9);
        }
        for (auto const& mesh : basic_mesh.meshes("cube"))
        {
            REQUIRE(mesh.topology() == element_topology::hexahedron8);
            REQUIRE(mesh.nodes_per_element() == 8);
            REQUIRE(mesh.elements() == 27);
        }
        for (auto const& mesh : basic_mesh.meshes("sides"))
        {
            REQUIRE(mesh.topology() == element_topology::quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 36);
        }
        for (auto const& mesh : basic_mesh.meshes("top"))
        {
            REQUIRE(mesh.topology() == element_topology::quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 9);
        }
    }
    SECTION("Test unique connectivities")
    {
        for (auto const& mesh : basic_mesh.meshes("bottom"))
        {
            auto const unique_node_list = mesh.unique_connectivity();

            std::vector<std::int32_t> const
                known_unique{0, 1, 2, 3, 14, 15, 16, 17, 18, 19, 26, 27, 36, 37, 38, 39};

            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("cube"))
        {
            auto const unique_node_list = mesh.unique_connectivity();

            REQUIRE(view::set_symmetric_difference(unique_node_list, view::ints(0, 64)).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("sides"))
        {
            auto const unique_node_list = mesh.unique_connectivity();

            std::vector<std::int32_t> const known_unique{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                         10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                         20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                         30, 31, 32, 33, 34, 35, 40, 41, 42, 43,
                                                         48, 49, 50, 51, 52, 53, 54, 55};

            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("top"))
        {
            auto const unique_node_list = mesh.unique_connectivity();

            std::vector<std::int32_t> const
                known_unique{4, 5, 6, 7, 10, 11, 22, 23, 28, 29, 30, 31, 44, 45, 46, 47};
            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }
    }
}
TEST_CASE("Solid submesh test")
{
    using namespace mechanical::solid;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh

    // Create the test objects
    basic_mesh basic_mesh(json::parse(cube_mesh_json()));
    nodal_coordinates nodal_coordinates(json::parse(cube_mesh_json()));

    auto& submeshes = basic_mesh.meshes("cube");

    REQUIRE(submeshes.size() == 1);

    auto& submesh = submeshes[0];

    auto mesh_coordinates = std::make_shared<material_coordinates>(nodal_coordinates.coordinates());

    mechanical::solid::fem_submesh fem_submesh(json::parse(material_data_json()),
                                               json::parse(simulation_data_json()),
                                               mesh_coordinates,
                                               submesh);

    int constexpr number_of_nodes = 64;
    int constexpr number_of_dofs = number_of_nodes * 3;
    int constexpr number_of_local_dofs = 8 * 3;

    vector displacement = 0.001 * vector::Random(number_of_dofs);

    mesh_coordinates->update_current_configuration(displacement);

    fem_submesh.update_internal_variables();

    auto& internal_vars = fem_submesh.internal_variables();

    SECTION("Degree of freedom test")
    {
        // Check the shape functions and dof definitions are ok
        REQUIRE(fem_submesh.dofs_per_node() == 3);
        REQUIRE(fem_submesh.shape_function().nodes() == 8);
    }
    SECTION("Default internal variables test")
    {
        // Check the standard ones are used
        REQUIRE(internal_vars.has(
            mechanical::solid::internal_variables_t::Tensor::DisplacementGradient));
        REQUIRE(
            internal_vars.has(mechanical::solid::internal_variables_t::Tensor::DeformationGradient));
        REQUIRE(internal_vars.has(mechanical::solid::internal_variables_t::Tensor::Cauchy));
        REQUIRE(internal_vars.has(mechanical::solid::internal_variables_t::Scalar::DetF));
    }
    SECTION("Tangent stiffness")
    {
        auto [local_dofs, stiffness] = fem_submesh.tangent_stiffness(0);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
        REQUIRE(stiffness.rows() == number_of_local_dofs);
        REQUIRE(stiffness.cols() == number_of_local_dofs);
        REQUIRE(stiffness.norm() != Approx(0.0).margin(ZERO_MARGIN));

        // Check symmetry for NeoHooke material model
        REQUIRE((stiffness - stiffness.transpose()).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Internal force")
    {
        auto [local_dofs, internal_force] = fem_submesh.internal_force(0);
        REQUIRE(internal_force.rows() == number_of_local_dofs);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
    }
    SECTION("Consistent and diagonal mass")
    {
        auto const& [local_dofs_0, mass_c] = fem_submesh.consistent_mass(0);
        auto const& [local_dofs_1, mass_d] = fem_submesh.diagonal_mass(0);

        REQUIRE(local_dofs_0.size() == number_of_local_dofs);
        REQUIRE(local_dofs_1.size() == number_of_local_dofs);

        REQUIRE(mass_c.rows() == number_of_local_dofs);
        REQUIRE(mass_c.cols() == number_of_local_dofs);

        REQUIRE(mass_d.rows() == number_of_local_dofs);
        REQUIRE(mass_d.cols() == 1);

        // Check the row sum is the same for each method
        for (auto i = 0; i < mass_d.rows(); i++)
        {
            REQUIRE(mass_c.row(i).sum() == Approx(mass_d(i)));
        }
    }
}
TEST_CASE("Solid mesh test")
{
    using mechanical::solid::fem_mesh;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    auto material_data = json::parse(material_data_json());
    auto simulation_data = json::parse(simulation_data_json());

    // Create the test objects
    basic_mesh basic_mesh(json::parse(cube_mesh_json()));
    nodal_coordinates nodal_coordinates(json::parse(cube_mesh_json()));

    REQUIRE(!simulation_data["Name"].empty());

    fem_mesh fem_mesh(basic_mesh,
                      material_data,
                      simulation_data,
                      simulation_data["Time"]["Increments"]["Initial"]);

    REQUIRE(fem_mesh.active_dofs() == 192);

    int constexpr number_of_nodes = 64;
    int constexpr number_of_dofs = number_of_nodes * 3;
    int constexpr number_of_local_dofs = 8 * 3;

    // Check that we only have one mesh group as we only have homogenous
    // element types
    REQUIRE(fem_mesh.meshes().size() == 1);

    for (auto const& fem_submesh : fem_mesh.meshes())
    {
        auto [local_dofs, internal_force] = fem_submesh.internal_force(0);
        REQUIRE(internal_force.rows() == number_of_local_dofs);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
    }

    SECTION("Check Dirichlet boundaries")
    {
        auto const& map = fem_mesh.displacement_boundaries();

        // See if correctly input in the map
        REQUIRE(map.find("bottom") != map.end());
        REQUIRE(map.find("top") != map.end());

        // And the negative is true...
        REQUIRE(map.find("sides") == map.end());
        REQUIRE(map.find("cube") == map.end());

        // Check the correct values in the boundary conditions
        for (auto const& fixed_bottom : map.find("bottom")->second)
        {
            REQUIRE(fixed_bottom.value_view(1.0) == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(fixed_bottom.dof_view().size() == 16);
        }

        for (auto const& disp_driven : map.find("top")->second)
        {
            REQUIRE(disp_driven.value_view(1.0) == Approx(0.001));
            REQUIRE(disp_driven.dof_view().size() == 16);
        }
    }
}
