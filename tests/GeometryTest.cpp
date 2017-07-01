
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

using namespace neon;

std::string cube_mesh_json()
{
    return "{\"Elements\":[{\"Indices\":[9,10,11,12,13,14,15,16,17],\"Name\":\"bottom\","
           "\"NodalConnectivity\":[[1,18,36,27],[27,36,37,26],[26,37,17,0],[18,19,38,36],[36,38,39,"
           "37],[37,39,16,17],[19,3,14,38],[38,14,15,39],[39,15,2,16]],\"Type\":3},{\"Indices\":["
           "67,80,79,78,77,76,75,74,73,72,71,70,69,68,54,66,65,64,63,62,61,60,59,58,57,56,55],"
           "\"Name\":\"cube\",\"NodalConnectivity\":[[62,60,56,58,63,61,57,59],[29,47,63,53,6,10,"
           "41,9],[28,45,62,52,29,47,63,53],[7,23,35,24,28,45,62,52],[47,46,61,63,10,11,40,41],[45,"
           "44,60,62,47,46,61,63],[23,22,34,35,45,44,60,62],[46,31,51,61,11,4,12,40],[44,30,50,60,"
           "46,31,51,61],[22,5,21,34,44,30,50,60],[53,63,59,55,9,41,43,8],[52,62,58,54,53,63,59,55]"
           ",[24,35,33,25,52,62,58,54],[63,61,57,59,41,40,42,43],[32,20,1,27,56,48,18,36],[35,34,"
           "32,33,62,60,56,58],[61,51,49,57,40,12,13,42],[60,50,48,56,61,51,49,57],[34,21,20,32,60,"
           "50,48,56],[55,59,39,16,8,43,15,2],[54,58,37,17,55,59,39,16],[25,33,26,0,54,58,37,17],["
           "59,57,38,39,43,42,14,15],[58,56,36,37,59,57,38,39],[33,32,27,26,58,56,36,37],[57,49,19,"
           "38,42,13,3,14],[56,48,18,36,57,49,19,38]],\"Type\":5},{\"Indices\":[45,37,38,39,40,41,"
           "42,43,44,36,46,47,48,49,50,51,52,53,18,1,2,3,4,5,6,7,8,0,19,20,21,22,23,24,25,26],"
           "\"Name\":\"sides\",\"NodalConnectivity\":[[7,24,52,28],[18,48,49,19],[19,49,13,3],[20,"
           "21,50,48],[48,50,51,49],[49,51,12,13],[21,5,30,50],[50,30,31,51],[51,31,4,12],[1,20,48,"
           "18],[28,52,53,29],[29,53,9,6],[24,25,54,52],[52,54,55,53],[53,55,8,9],[25,0,17,54],[54,"
           "17,16,55],[55,16,2,8],[4,12,40,11],[27,32,33,26],[26,33,25,0],[20,21,34,32],[32,34,35,"
           "33],[33,35,24,25],[21,5,22,34],[34,22,23,35],[35,23,7,24],[1,20,32,27],[11,40,41,10],["
           "10,41,9,6],[12,13,42,40],[40,42,43,41],[41,43,8,9],[13,3,14,42],[42,14,15,43],[43,15,2,"
           "8]],\"Type\":3},{\"Indices\":[27,28,29,30,31,32,33,34,35],\"Name\":\"top\","
           "\"NodalConnectivity\":[[5,30,44,22],[22,44,45,23],[23,45,28,7],[30,31,46,44],[44,46,47,"
           "45],[45,47,29,28],[31,4,11,46],[46,11,10,47],[47,10,6,29]],\"Type\":3}],\"Nodes\":[{"
           "\"Coordinates\":[[0,0,0],[1,0,0],[0,1,0],[1,1,0],[1,1,1],[1,0,1],[0,1,1],[0,0,1],[0,1,"
           "0.33333333333250098],[0,1,0.66666666666578744],[0.33333333333250098,1,1],[0."
           "66666666666578744,1,1],[1,1,0.66666666666759178],[1,1,0.33333333333472071],[0."
           "66666666666759178,1,0],[0.33333333333472071,1,0],[0,0.66666666666759178,0],[0,0."
           "33333333333472071,0],[1,0.33333333333250098,0],[1,0.66666666666578744,0],[1,0,0."
           "33333333333250098],[1,0,0.66666666666578744],[0.66666666666759178,0,1],[0."
           "33333333333472071,0,1],[0,0,0.66666666666759178],[0,0,0.33333333333472071],[0."
           "33333333333250098,0,0],[0.66666666666578744,0,0],[0,0.33333333333250098,1],[0,0."
           "66666666666578744,1],[1,0.33333333333250098,1],[1,0.66666666666578744,1],[0."
           "66666666666638896,0,0.33333333333324089],[0.333333333333241,0,0.33333333333398069],[0."
           "66666666666699048,0,0.66666666666638896],[0.33333333333398091,0,0.66666666666699026],["
           "0.66666666666638896,0.33333333333324089,0],[0.333333333333241,0.33333333333398069,0],["
           "0.66666666666699048,0.66666666666638896,0],[0.33333333333398091,0.66666666666699026,0],"
           "[0.66666666666638896,1,0.66666666666699026],[0.333333333333241,1,0.66666666666638896],["
           "0.66666666666699048,1,0.33333333333398069],[0.33333333333398091,1,0.333333333333241],["
           "0.66666666666699026,0.33333333333250098,1],[0.33333333333398091,0.33333333333250098,1],"
           "[0.66666666666638874,0.66666666666578744,1],[0.33333333333324089,0.66666666666578744,1]"
           ",[1,0.33333333333250098,0.33333333333324089],[1,0.66666666666578744,0."
           "33333333333398091],[1,0.33333333333250098,0.66666666666638874],[1,0.66666666666578744,"
           "0.66666666666699048],[0,0.333333333333241,0.66666666666699026],[0,0.66666666666638896,"
           "0.66666666666638896],[0,0.33333333333398091,0.33333333333398091],[0,0."
           "66666666666699048,0.333333333333241],[0.66666666666658936,0.33333333333299442,0."
           "33333333333348758],[0.66666666666679042,0.66666666666618857,0.33333333333373422],[0."
           "33333333333348752,0.33333333333348752,0.33333333333373438],[0.33333333333373433,0."
           "66666666666658925,0.33333333333348758],[0.66666666666678998,0.33333333333274773,0."
           "66666666666658947],[0.66666666666658991,0.66666666666598795,0.66666666666679009],[0."
           "33333333333373422,0.33333333333299431,0.66666666666678975],[0.33333333333348769,0."
           "66666666666618846,0.66666666666658936]],\"Indices\":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,"
           "14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,"
           "43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]}]}";
}

std::string material_data_json()
{
    return "{\"Name\" : \"steel\",\"ElasticModulus\" : 200.0e9, \"PoissonsRatio\" "
           ": 0.3 }";
}

std::string simulation_data_json()
{
    return "{ \"BoundaryConditions\" : [ {\"Name\" : \"bottom\", "
           "\"Type\" : "
           "\"Displacement\", \"Values\" : {\"x\" : 0.0, \"y\" : 0.0, \"z\" : 0.0}}, "
           "{\"Name\" : "
           "\"top\", \"Type\" : \"Displacement\", \"Values\" : {\"z\" : 0.001}} ], "
           "\"ConstitutiveModel\" : \"NeoHooke\", "
           "\"ElementOptions\" : {\"Quadrature\" : \"Reduced\"}, "
           "\"Name\" : \"cube\"}";
}

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
    List local_node_list = ranges::view::ints(0, 3);

    // Check that we are fetching the right local element vector
    List local_dof_list = ranges::view::ints(0, 9);

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
        REQUIRE((material_coordinates.displacement() - local_displacements).norm() == Approx(0.0));
    }
    SECTION("Test local displacement via lookup")
    {
        material_coordinates.update_current_configuration(local_displacements);
        REQUIRE((material_coordinates.displacement(local_dof_list) - local_displacements).norm() ==
                Approx(0.0));
    }
    SECTION("Test element view configuration")
    {
        REQUIRE((material_coordinates.initial_configuration(local_node_list) - local_initial_config)
                    .norm() == Approx(0.0));
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

            REQUIRE(ranges::view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("cube"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            REQUIRE(ranges::view::set_symmetric_difference(unique_node_list, ranges::view::ints(0, 64))
                        .empty());
        }

        for (auto const& mesh : basic_mesh.meshes("sides"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            List const known_unique{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
                                    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                                    32, 33, 34, 35, 40, 41, 42, 43, 48, 49, 50, 51, 52, 53, 54, 55};

            REQUIRE(ranges::view::set_symmetric_difference(unique_node_list, known_unique).empty());
        }

        for (auto const& mesh : basic_mesh.meshes("top"))
        {
            auto const unique_node_list = mesh.unique_connectivities();

            List const known_unique{4, 5, 6, 7, 10, 11, 22, 23, 28, 29, 30, 31, 44, 45, 46, 47};
            REQUIRE(ranges::view::set_symmetric_difference(unique_node_list, known_unique).empty());
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

    auto material_coordinates =
        std::make_shared<MaterialCoordinates>(nodal_coordinates.coordinates());

    solid::femSubmesh fem_submesh(material_data, simulation_data, material_coordinates, submesh);

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
        REQUIRE(internal_vars.has(InternalVariables::Tensor::CauchyStress));
        REQUIRE(internal_vars.has(InternalVariables::Scalar::DetF));
    }
    SECTION("Material model internal variables")
    {
        REQUIRE(internal_vars.has(InternalVariables::Tensor::Kirchhoff));
    }
    SECTION("Tangent stiffness")
    {
        auto[local_dofs, stiffness] = fem_submesh.tangent_stiffness(0);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
        REQUIRE(stiffness.rows() == number_of_local_dofs);
        REQUIRE(stiffness.cols() == number_of_local_dofs);
        REQUIRE(stiffness.norm() != Approx(0.0));
    }
    SECTION("Internal force")
    {
        auto[local_dofs, internal_force] = fem_submesh.internal_force(0);
        REQUIRE(internal_force.rows() == number_of_local_dofs);
        REQUIRE(local_dofs.size() == number_of_local_dofs);
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
        auto const& map = fem_mesh.dirichlet_boundary_map();

        // See if correctly input in the map
        REQUIRE(map.find("bottom") != map.end());
        REQUIRE(map.find("top") != map.end());

        // And the negative is true...
        REQUIRE(map.find("sides") == map.end());
        REQUIRE(map.find("cube") == map.end());

        // Check the correct values in the boundary conditions
        for (auto const& fixed_bottom : map.find("bottom")->second)
        {
            REQUIRE(fixed_bottom.prescribed_value() == Approx(0.0));
            REQUIRE(fixed_bottom.prescribed_dofs().size() == 16);
        }

        for (auto const& disp_driven : map.find("top")->second)
        {
            REQUIRE(disp_driven.prescribed_value() == Approx(0.001));
            REQUIRE(disp_driven.prescribed_dofs().size() == 16);
        }
    }
}
TEST_CASE("Nodal ordering Adapater")
{
    NodeOrderingAdapter node_adapter;

    SECTION("Tetrahedron10 ordering")
    {
        // Create a dummy connectivity set
        std::vector<List> nodal_connectivity = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                                {10, 11, 12, 13, 14, 15, 16, 17, 18, 19}};

        node_adapter.convert_from_gmsh(nodal_connectivity, ElementTopology::Tetrahedron10);

        REQUIRE(nodal_connectivity[0][4] == 9);
        REQUIRE(nodal_connectivity[1][4] == 19);

        REQUIRE(nodal_connectivity[0][9] == 4);
        REQUIRE(nodal_connectivity[1][9] == 14);
    }
}
