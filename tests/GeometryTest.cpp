
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"
#include "mesh/ElementTopology.hpp"
#include "mesh/NodeOrderingAdapter.hpp"

#include "mesh/solid/MaterialCoordinates.hpp"
#include "mesh/solid/Submesh.hpp"
#include "mesh/solid/femMesh.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;
using namespace ranges;

TEST_CASE("Testing material coordinates", "[MaterialCoordinates]")
{
    // Build a right angled triangle
    Vector coordinates(9);
    coordinates << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    // Setup the test case
    NodalCoordinates nodes(coordinates);
    solid::MaterialCoordinates material_coordinates(coordinates);

    // Test with a random displacement vector
    Vector local_displacements = Vector::Random(9);

    Matrix local_initial_config(3, 3);
    local_initial_config << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    // Indices for the first three nodes
    List local_node_list = view::ints(0, 3);

    // Check that we are fetching the right local element vector
    List local_dof_list = view::ints(0, 9);

    SECTION("Nodes scaffolding")
    {
        Vector triangle = Vector::Zero(9);
        triangle << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
        REQUIRE((nodes.coordinates() - triangle).norm() == Approx(0.0));
    }
    SECTION("Check initial displacements are zero")
    {
        REQUIRE(material_coordinates.displacement().norm() == Approx(0.0));
    }
    SECTION("Test update of coordinates")
    {
        material_coordinates.update_current_configuration(local_displacements);
        REQUIRE((material_coordinates.displacement() - local_displacements).norm()
                == Approx(0.0));
    }
    SECTION("Test local displacement via lookup")
    {
        material_coordinates.update_current_configuration(local_displacements);
        REQUIRE(
            (material_coordinates.displacement(local_dof_list) - local_displacements).norm()
            == Approx(0.0));
    }
    SECTION("Test element view configuration")
    {
        REQUIRE((material_coordinates.initial_configuration(local_node_list)
                 - local_initial_config)
                    .norm()
                == Approx(0.0));
    }
}
TEST_CASE("Basic mesh test")
{
    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data;
    Json::Reader mesh_file;

    REQUIRE(mesh_file.parse(cube_mesh_json().c_str(), mesh_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);

    int constexpr number_of_nodes = 64;

    REQUIRE(nodal_coordinates.size() == number_of_nodes);

    SECTION("Test corner vertices")
    {
        REQUIRE((nodal_coordinates[{0}] - Vector3(0.0, 0.0, 0.0)).norm() == Approx(0.0));
        REQUIRE((nodal_coordinates[{1}] - Vector3(1.0, 0.0, 0.0)).norm() == Approx(0.0));
        REQUIRE((nodal_coordinates[{2}] - Vector3(0.0, 1.0, 0.0)).norm() == Approx(0.0));
        REQUIRE((nodal_coordinates[{3}] - Vector3(1.0, 1.0, 0.0)).norm() == Approx(0.0));
    }
    SECTION("Test mesh data for boundary and volume elements")
    {
        for (auto const& mesh : basic_mesh.meshes("bottom"))
        {
            REQUIRE(mesh.topology() == ElementTopology::Quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 9);
        }
        for (auto const& mesh : basic_mesh.meshes("cube"))
        {
            REQUIRE(mesh.topology() == ElementTopology::Hexahedron8);
            REQUIRE(mesh.nodes_per_element() == 8);
            REQUIRE(mesh.elements() == 27);
        }
        for (auto const& mesh : basic_mesh.meshes("sides"))
        {
            REQUIRE(mesh.topology() == ElementTopology::Quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 36);
        }
        for (auto const& mesh : basic_mesh.meshes("top"))
        {
            REQUIRE(mesh.topology() == ElementTopology::Quadrilateral4);
            REQUIRE(mesh.nodes_per_element() == 4);
            REQUIRE(mesh.elements() == 9);
        }
    }
    SECTION("Test unique connectivities")
    {
        for (auto const& mesh : basic_mesh.meshes("bottom"))
        {
            auto const unique_node_list = mesh.unique_connectivities();
            List const known_unique{0, 1, 2, 3, 14, 15, 16, 17, 18, 19, 26, 27, 36, 37, 38, 39};

            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("cube"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            REQUIRE(
                view::set_symmetric_difference(unique_node_list, view::ints(0, 64)).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("sides"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            List const known_unique{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                                    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                    24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                    40, 41, 42, 43, 48, 49, 50, 51, 52, 53, 54, 55};

            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("top"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            List const known_unique{4, 5, 6, 7, 10, 11, 22, 23, 28, 29, 30, 31, 44, 45, 46, 47};
            REQUIRE(view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }
    }
}
TEST_CASE("Solid submesh test")
{
    using solid::MaterialCoordinates;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data;
    Json::Reader mesh_file, material_file, simulation_file;

    REQUIRE(mesh_file.parse(cube_mesh_json().c_str(), mesh_data));
    REQUIRE(material_file.parse(material_data_json().c_str(), material_data));
    REQUIRE(simulation_file.parse(simulation_data_json().c_str(), simulation_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);

    auto& submeshes = basic_mesh.meshes("cube");

    REQUIRE(submeshes.size() == 1);

    auto& submesh = submeshes[0];

    auto material_coordinates = std::make_shared<MaterialCoordinates>(
        nodal_coordinates.coordinates());

    solid::femSubmesh fem_submesh(material_data,
                                  simulation_data,
                                  material_coordinates,
                                  submesh);

    int constexpr number_of_nodes = 64;
    int constexpr number_of_dofs = number_of_nodes * 3;
    int constexpr number_of_local_dofs = 8 * 3;

    Vector displacement = 0.001 * Vector::Random(number_of_dofs);

    material_coordinates->update_current_configuration(displacement);

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
        REQUIRE(internal_vars.has(InternalVariables::Tensor::DisplacementGradient));
        REQUIRE(internal_vars.has(InternalVariables::Tensor::DeformationGradient));
        REQUIRE(internal_vars.has(InternalVariables::Tensor::Cauchy));
        REQUIRE(internal_vars.has(InternalVariables::Scalar::DetF));
    }
    SECTION("Tangent stiffness")
    {
        auto[local_dofs, stiffness] = fem_submesh.tangent_stiffness(0);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
        REQUIRE(stiffness.rows() == number_of_local_dofs);
        REQUIRE(stiffness.cols() == number_of_local_dofs);
        REQUIRE(stiffness.norm() != Approx(0.0));

        // Check symmetry for NeoHooke material model
        REQUIRE((stiffness - stiffness.transpose()).norm() == Approx(0.0));
    }
    SECTION("Internal force")
    {
        auto[local_dofs, internal_force] = fem_submesh.internal_force(0);
        REQUIRE(internal_force.rows() == number_of_local_dofs);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
    }
    SECTION("Consistent and diagonal mass")
    {
        auto const & [ local_dofs_0, mass_c ] = fem_submesh.consistent_mass(0);
        auto const & [ local_dofs_1, mass_d ] = fem_submesh.diagonal_mass(0);

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
    using solid::MaterialCoordinates;
    using solid::femMesh;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data;
    Json::Reader mesh_file, material_file, simulation_file;

    REQUIRE(mesh_file.parse(cube_mesh_json().c_str(), mesh_data));
    REQUIRE(material_file.parse(material_data_json().c_str(), material_data));
    REQUIRE(simulation_file.parse(simulation_data_json().c_str(), simulation_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);

    REQUIRE(!simulation_data["Name"].empty());

    femMesh fem_mesh(basic_mesh, material_data, simulation_data);

    REQUIRE(fem_mesh.active_dofs() == 192);

    int constexpr number_of_nodes = 64;
    int constexpr number_of_dofs = number_of_nodes * 3;
    int constexpr number_of_local_dofs = 8 * 3;

    // Check that we only have one mesh group as we only have homogenous
    // element types
    REQUIRE(fem_mesh.meshes().size() == 1);

    for (auto const& fem_submesh : fem_mesh.meshes())
    {
        auto[local_dofs, internal_force] = fem_submesh.internal_force(0);
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
            REQUIRE(fixed_bottom.value_view(1.0) == Approx(0.0));
            REQUIRE(fixed_bottom.dof_view().size() == 16);
        }

        for (auto const& disp_driven : map.find("top")->second)
        {
            REQUIRE(disp_driven.value_view(1.0) == Approx(0.001));
            REQUIRE(disp_driven.dof_view().size() == 16);
        }
    }
}
TEST_CASE("Nodal ordering Adapater")
{
    NodeOrderingAdapter node_adapter;

    SECTION("Tetrahedron10 ordering")
    {
        // Create a dummy connectivity set
        std::vector<List> nodal_connectivity = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};

        node_adapter.convert_from_gmsh(nodal_connectivity, ElementTopology::Tetrahedron10);

        REQUIRE(nodal_connectivity[0][4] == 9);
        REQUIRE(nodal_connectivity[0][9] == 4);
        REQUIRE(nodal_connectivity[0][3] == 0);
        REQUIRE(nodal_connectivity[0][0] == 3);
    }
}
