
#include "SimulationControl.hpp"

#include "Exceptions.hpp"
#include "mesh/solid/femMesh.hpp"

#include "assembler/solid/femStaticMatrix.hpp"

#include "assembler/solid/femDynamicMatrix.hpp"

#include <iomanip>
#include <thread>

#include <boost/filesystem.hpp>
#include <json/reader.h>
#include <termcolor/termcolor.hpp>

namespace neon
{
int SimulationControl::threads = std::thread::hardware_concurrency();

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
    auto start = std::chrono::high_resolution_clock::now();

    this->print_banner();

    std::cout << std::string(2, ' ') << termcolor::bold << "Preprocessing mesh and simulation data\n"
              << termcolor::reset;

    std::ifstream file(input_file_name);

    Json::Reader reader;

    if (!reader.parse(file, root, false))
    {
        throw std::runtime_error(reader.getFormattedErrorMessages());
    }

    this->check_input_fields(root);

    if (root.isMember("Cores")) threads = root["Cores"].asInt();

    auto const material_names = this->parse_material_names(root["Material"]);
    auto const part_names = this->parse_part_names(root["Part"], material_names);

    // Add in the parts and populate the mesh stores
    for (auto const& part : root["Part"])
    {
        Json::Value mesh_file;
        Json::Reader mesh_reader;

        auto const material = *ranges::find_if(root["Material"], [&part](auto const& material) {
            return material["Name"].asString() == part["Material"].asString();
        });

        std::ifstream mesh_input_stream(part["Name"].asString() + ".mesh");

        if (!mesh_reader.parse(mesh_input_stream, mesh_file, false))
        {
            throw std::runtime_error(mesh_reader.getFormattedErrorMessages());
        }

        mesh_store.try_emplace(part["Name"].asString(), mesh_file, material);

        std::cout << std::string(4, ' ') << "Inserted " << part["Name"].asString()
                  << " into the mesh store\n";
    }

    std::array<std::string, 6> const required_fields = {
        {"Name", "Time", "Solution", "Visualisation", "LinearSolver", "NonlinearOptions"}};

    // Build a list of all the load steps for a given mesh
    for (auto const& simulation : root["SimulationCases"])
    {
        // Ensure the required fields exist
        for (auto const& required_field : required_fields)
        {
            if (!simulation.isMember(required_field))
                throw std::runtime_error("A simulation case needs a \"" + required_field
                                         + "\" field\n");
        }

        // Multibody simulations not (yet) supported
        assert(simulation["Mesh"].size() == 1);

        // Make sure the simulation mesh exists in the mesh store
        if (mesh_store.find(simulation["Mesh"][0]["Name"].asString()) == mesh_store.end())
        {
            throw std::runtime_error("Mesh name \"" + simulation["Mesh"][0]["Name"].asString()
                                     + "\" was not found in the mesh store");
        }
    }

    build_simulation_tree();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << termcolor::bold << termcolor::green << std::string(2, ' ')
              << "Preprocessing complete in " << elapsed_seconds.count() << "s\n"
              << termcolor::reset << std::endl;
}

void SimulationControl::start()
{
    // It is desirable to have multiple load steps throughout the simulation.
    // This can be achieved by stepping through each simulation and performing
    // an internal restart, using the data from the previous time step.
    for (auto const& simulation : root["SimulationCases"])
    {
        if (!simulation["Inherits"].empty()) continue;

        std::cout << std::string(2, ' ') << termcolor::bold << "Simulating case \""
                  << simulation["Name"].asString() << "\"" << termcolor::reset << std::endl
                  << std::endl;

        auto const& mesh_data = simulation["Mesh"][0];

        auto simulation_mesh = mesh_store.find(mesh_data["Name"].asString());

        std::cout << std::string(4, ' ') << "Module \"" << simulation["Type"].asString() << "\"\n";

        std::cout << std::string(4, ' ') << "Solution \"" << simulation["Solution"].asString()
                  << "\"\n";

        auto const & [ mesh, material ] = simulation_mesh->second;

        solid::femMesh fem_mesh(mesh, material, mesh_data);

        solid::femStaticMatrix fem_matrix(fem_mesh,
                                          Visualisation(root["Name"].asString(),
                                                        fem_mesh,
                                                        simulation["Visualisation"]),
                                          simulation["LinearSolver"],
                                          simulation["NonlinearOptions"],
                                          simulation["Time"]);

        fem_matrix.solve();

        for (auto const& next_step : multistep_simulations[simulation["Name"].asString()])
        {
            std::cout << termcolor::bold << "\n"
                      << std::string(2, ' ') << "Simulating case \"" << next_step["Name"].asString()
                      << "\"" << termcolor::reset << std::endl
                      << std::endl;

            fem_mesh.internal_restart(next_step["Mesh"][0]);

            fem_matrix.internal_restart(next_step["LinearSolver"], next_step["Time"]);

            fem_matrix.solve();
        }
    }
}

void SimulationControl::build_simulation_tree()
{
    // For each simulation step, the total number of subsequent relationships
    // need to be determined, such that an analysis can be performed in order
    // Build a list of all the load steps for a given mesh
    for (auto const& simulation : root["SimulationCases"])
    {
        if (simulation["Inherits"].empty())
        {
            find_children(simulation["Name"].asString(), simulation["Name"].asString());
        }
    }

    for (auto const & [ name, queue ] : multistep_simulations)
    {
        std::cout << std::string(4, ' ') << "Simulation \"" << name << "\" is continued by:\n";
        for (auto const& item : queue)
        {
            std::cout << "\t\"" << item["Name"].asString() << "\"" << std::endl;
        }
    }
}

void SimulationControl::find_children(std::string const& parent_name,
                                      std::string const& next_parent_name)
{
    for (auto const& simulation : root["SimulationCases"])
    {
        if (simulation["Inherits"].empty()) continue;

        if (simulation["Inherits"].asString() == next_parent_name)
        {
            multistep_simulations[parent_name].push_back(simulation);
            // Recursive step
            find_children(parent_name, simulation["Name"].asString());
        }
    }
}

void SimulationControl::print_banner() const
{
    std::string const welcome_message("neon - a non-linear finite element code");

    std::cout << termcolor::bold;
    std::cout << std::setw(welcome_message.length() + 9) << std::setfill('=') << "\n";
    std::cout << std::string(4, ' ') << welcome_message << "\n";
    std::cout << std::setw(welcome_message.length() + 9) << std::setfill('=') << "\n";
    std::cout << termcolor::reset << std::endl << std::setfill(' ');
}

void SimulationControl::check_input_fields(Json::Value const& root) const
{
    // Check the important fields exist before anything else is done
    if (!root.isMember("Part")) throw std::runtime_error("Part is not in input file");
    if (!root.isMember("Name")) throw std::runtime_error("Name is not in input file");
    if (!root.isMember("Material")) throw std::runtime_error("Material is not in input file");
    if (!root.isMember("SimulationCases"))
        throw std::runtime_error("SimulationCases is not in input file");
}

/** Extract the material names from the input file */
std::unordered_set<std::string> SimulationControl::parse_material_names(Json::Value const& materials) const
{
    std::unordered_set<std::string> material_names;

    for (auto const& material : root["Material"])
    {
        if (material["Name"].empty()) throw std::runtime_error("Material: Name");

        auto const[it, inserted] = material_names.emplace(material["Name"].asString());

        if (!inserted) throw DuplicateNameException("Material");
    }
    return material_names;
}

/** Extract the part names from the input file */
std::unordered_set<std::string> SimulationControl::parse_part_names(
    Json::Value const& parts, std::unordered_set<std::string> const& material_names) const
{
    std::unordered_set<std::string> part_names;

    // Load in all the part names and error check
    for (auto const& part : root["Part"])
    {
        if (part["Name"].empty()) throw std::runtime_error("Part: Name");

        if (ranges::find(material_names, part["Material"].asString()) == material_names.end())
        {
            throw std::runtime_error("The part material was not found in the provided "
                                     "materials\n");
        }

        auto const[it, inserted] = part_names.emplace(part["Name"].asString());

        if (!inserted) throw DuplicateNameException("Part");
    }
    return part_names;
}
}
