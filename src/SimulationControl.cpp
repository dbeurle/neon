
#include "SimulationControl.hpp"

#include "PreprocessorExceptions.hpp"

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

    std::cout << "=============================================================\n";
    std::cout << " neon - a non-linear finite element code\n\n";
    std::cout << "    Simulation to use " << threads << " thread(s)\n";
    std::cout << "    Preprocessing for gmsh input completed.\n";
    std::cout << "*------------------------------------------------------------*\n";

    std::cout << "Mesh name " << root["Part"][0]["Name"].asString() << std::endl;

    solid::femMesh fem_mesh(mesh_store.find(root["Part"][0]["Name"].asString())->second,
                            root["Material"][0],
                            root["Simulation"][0]["Mesh"][0]);

    solid::femMatrix fem_matrix(fem_mesh, root["Simulation"][0]["LinearSolver"]);

    fem_matrix.solve();
}

void SimulationControl::start() {}
}
