
#include "SimulationControl.hpp"

#include "PreprocessorExceptions.hpp"

#include "assembler/solid/femStaticMatrix.hpp"

#include "assembler/solid/femDynamicMatrix.hpp"

#include "mesh/solid/femMesh.hpp"

#include <unordered_set>

#include <boost/filesystem.hpp>
#include <json/json.h>

namespace neon
{
SimulationControl::SimulationControl(std::string const& input_file_name)
    : input_file_name(input_file_name)
{
    if (input_file_name == "" || input_file_name.empty()) throw NoInputException();

    boost::filesystem::path input_path(input_file_name);

    // Strip the extension from the filename
    std::string extension = boost::filesystem::extension(input_path);
    std::string base_name = boost::filesystem::basename(input_path);

    // Attempt to open the json input file
    if (extension != ".json") throw InvalidExtensionException(extension);

    this->parse();
}

void SimulationControl::parse()
{
    std::cout << "=============================================================\n\n";
    std::cout << "           neon - a non-linear finite element code\n\n";
    std::cout << "=============================================================\n\n";

    std::ifstream file(input_file_name);

    Json::Value root;
    Json::Reader reader;

    if (!reader.parse(file, root, false))
        throw JsonFileParseException(reader.getFormattedErrorMessages());

    // Check the important fields exist before anything else is done
    if (root["Part"].empty()) throw EmptyFieldException("Part");
    if (root["Name"].empty()) throw EmptyFieldException("Name");
    if (root["Material"].empty()) throw EmptyFieldException("Material");
    if (root["Simulation"].empty()) throw EmptyFieldException("Simulation");

    // auto const& simulation_name = root["AnalysisName"].asString();
    if (!root["Cores"].empty()) threads = root["Cores"].asInt();

    std::unordered_set<std::string> material_names, part_names;

    // Load in all the material names
    for (const auto& material : root["Material"])
    {
        if (material["Name"].empty()) throw EmptyFieldException("Material: Name");

        auto const[it, inserted] = material_names.emplace(material["Name"].asString());

        if (!inserted) throw DuplicateNameException("Material");
    }

    // Load in all the part names
    for (const auto& part : root["Part"])
    {
        if (part["Name"].empty()) throw EmptyFieldException("Part: Name");

        auto const[it, inserted] = part_names.emplace(part["Name"].asString());

        if (!inserted) throw DuplicateNameException("Part");
    }

    // Add in the parts and populate the mesh stores
    for (const auto& part : root["Part"])
    {
        Json::Value mesh_file;
        Json::Reader mesh_reader;

        std::ifstream mesh_input_stream(part["Name"].asString() + ".mesh");

        if (!mesh_reader.parse(mesh_input_stream, mesh_file, false))
        {
            throw JsonFileParseException(mesh_reader.getFormattedErrorMessages());
        }
        mesh_store.try_emplace(part["Name"].asString(), mesh_file);
    }

    std::cout << "    Simulation to use " << threads << " thread(s)\n";
    std::cout << "    Preprocessing complete\n";

    std::cout << "Mesh name " << root["Part"][0]["Name"].asString() << std::endl;

    solid::femMesh fem_mesh(mesh_store.find(root["Part"][0]["Name"].asString())->second,
                            root["Material"][0],
                            root["Simulation"][0]["Mesh"][0]);

    if (root["Simulation"][0]["Solution"].empty())
    {
        throw std::runtime_error(
            "Simulations need a solution (\"Transient\" or \"Equilibrium\") type\n");
    }

    if (root["Simulation"][0]["Solution"].asString() == "Transient")
    {
        if (root["Simulation"][0]["Time"].empty())
            throw std::runtime_error("Simulation needs a \"Time\" field specified");

        solid::femDynamicMatrix fem_matrix(fem_mesh,
                                           root["Simulation"][0]["LinearSolver"],
                                           root["Simulation"][0]["Time"]);
        fem_matrix.solve();
    }
    else if (root["Simulation"][0]["Solution"].asString() == "Equilibrium")
    {
        if (root["Simulation"][0]["Time"].empty())
        {
            throw std::runtime_error("Simulations need a \"Time\" field\n");
        }

        solid::femStaticMatrix fem_matrix(fem_mesh,
                                          root["Simulation"][0]["LinearSolver"],
                                          root["Simulation"][0]["Time"]);
        fem_matrix.solve();

        std::cout << "*-------------------------------*\n";
        std::cout << "  Setting up the next load step\n";
        std::cout << "*-------------------------------*\n";

        // Update the mesh with a continuation
        fem_mesh.continuation(root["Simulation"][1]["Mesh"][0]);

        // Update the matrix with the new time solution
        fem_matrix.continuation(root["Simulation"][1]["Time"]);

        fem_matrix.solve();
    }
    else
    {
        std::runtime_error("Specify solution type as Transient or Equilibrium\n");
    }
}

void SimulationControl::start() {}
}
