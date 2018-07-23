
#include "simulation_parser.hpp"

#include "exceptions.hpp"
#include "modules/abstract_module.hpp"
#include "modules/module_factory.hpp"

#include <iomanip>
#include <fstream>
#include <thread>

#include <boost/filesystem.hpp>
#include <termcolor/termcolor.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/find_if.hpp>

namespace neon
{
int simulation_parser::threads = std::thread::hardware_concurrency();

simulation_parser::simulation_parser(std::string const& input_file_name)
    : input_file_name(input_file_name)
{
    if (input_file_name == "" || input_file_name.empty())
    {
        throw std::domain_error("No input file found.  An input file needs to be provided: "
                                "\"<filename>.neon\"\n");
    }

    boost::filesystem::path input_path(input_file_name);

    // Strip the extension from the filename
    std::string extension = boost::filesystem::extension(input_path);
    std::string base_name = boost::filesystem::basename(input_path);

    if (extension != ".json")
    {
        throw std::domain_error("Extension \"" + extension + "\" is not supported, use \".json\"");
    }
    this->parse();
}

simulation_parser::~simulation_parser() = default;

void simulation_parser::parse()
{
    auto const start = std::chrono::steady_clock::now();

    this->print_banner();

    std::cout << std::string(2, ' ') << termcolor::bold << "Preprocessing mesh and simulation data\n"
              << termcolor::reset;

    {
        std::ifstream file(input_file_name);
        file >> root;
    }

    this->check_input_fields();

    if (root.find("Cores") != root.end())
    {
        threads = root["Cores"];
    }

    auto const material_names = this->parse_material_names(root["Material"]);
    auto const part_names = this->parse_part_names(root["Part"], material_names);

    // Add in the parts and populate the mesh stores
    for (auto const& part : root["Part"])
    {
        auto const material = *ranges::find_if(root["Material"], [&part](auto const& material) {
            return material["Name"] == part["Material"];
        });

        auto const read_start = std::chrono::steady_clock::now();

        std::ifstream mesh_input_stream(part["Name"].get<std::string>() + ".mesh");

        json mesh_file;
        mesh_input_stream >> mesh_file;

        auto const read_end = std::chrono::steady_clock::now();

        std::cout << std::string(4, ' ') << "Parsed " << part["Name"] << " mesh from file in "
                  << std::chrono::duration<double>(read_end - read_start).count() << "s\n";

        mesh_store.try_emplace(part["Name"], mesh_file, material);

        std::cout << std::string(4, ' ') << "Allocated internal storage for " << part["Name"]
                  << " in "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - read_end).count()
                  << "s\n";
    }

    // Build a list of all the load steps for a given mesh
    for (auto const& simulation : root["SimulationCases"])
    {
        // Ensure the required fields exist
        for (auto required_field : {"Name", "Time", "Solution", "LinearSolver"})
        {
            if (!simulation.count(required_field))
            {
                throw std::domain_error("A simulation case needs a \"" + std::string(required_field)
                                        + "\" field\n");
            }
        }
        // Multibody simulations not (yet) supported
        assert(simulation["Mesh"].size() == 1);

        // Make sure the simulation mesh exists in the mesh store
        if (mesh_store.find(simulation["Mesh"][0]["Name"]) == mesh_store.end())
        {
            throw std::domain_error("Mesh name \"" + simulation["Mesh"][0]["Name"].get<std::string>()
                                    + "\" was not found in the mesh store");
        }
    }
    build_simulation_tree();

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << termcolor::bold << termcolor::green << std::string(2, ' ')
              << "Preprocessing complete in " << elapsed_seconds.count() << "s\n"
              << termcolor::reset << std::endl;
}

void simulation_parser::start()
{
    // Construct the modules and automatically check for correct input.
    // Throws the appropriate exception when an error is detected.
    for (auto const& [name, simulations] : multistep_simulations)
    {
        for (auto const& simulation : simulations)
        {
            modules.emplace_back(make_module(simulation, mesh_store));
        }
    }
    for (auto const& module : modules)
    {
        module->perform_simulation();
    }
}

void simulation_parser::build_simulation_tree()
{
    // For each simulation step, the total number of subsequent relationships
    // need to be determined, such that an analysis can be performed in order
    // Build a list of all the load steps for a given mesh
    for (auto const& simulation : root["SimulationCases"])
    {
        if (simulation.find("Inherits") == simulation.end())
        {
            std::string const& simulation_name = simulation["Name"];

            multistep_simulations[simulation_name].emplace_front(simulation);

            find_children(simulation_name, simulation_name);
        }
    }

    // Print out continuation data if applicable
    for (auto const& [name, queue] : multistep_simulations)
    {
        if (!queue.empty())
        {
            std::cout << std::string(4, ' ') << "Simulation \"" << name << "\" is continued by:\n";

            for (auto const& item : queue)
            {
                std::cout << std::string(4, ' ') << item["Name"] << std::endl;
            }
        }
    }
}

void simulation_parser::find_children(std::string const& parent_name,
                                      std::string const& next_parent_name)
{
    for (auto const& simulation : root["SimulationCases"])
    {
        if (simulation.find("Inherits") != simulation.end()) continue;

        if (simulation["Inherits"].get<std::string>() == next_parent_name)
        {
            multistep_simulations[parent_name].push_back(simulation);
            // Recursive step
            find_children(parent_name, simulation["Name"]);
        }
    }
}

void simulation_parser::print_banner() const
{
    std::string const welcome_message("neon - a non-linear finite element code");

    std::cout << termcolor::bold;
    std::cout << std::setw(welcome_message.length() + 9) << std::setfill('=') << "\n";
    std::cout << std::string(4, ' ') << welcome_message << "\n";
    std::cout << std::setw(welcome_message.length() + 9) << std::setfill('=') << "\n";
    std::cout << termcolor::reset << std::endl << std::setfill(' ');
}

void simulation_parser::check_input_fields() const
{
    // Check the important fields exist before anything else is done
    if (root.find("Part") == root.end())
    {
        throw std::domain_error("Part is not in input file");
    }
    if (root.find("Name") == root.end())
    {
        throw std::domain_error("Name is not in input file");
    }
    if (root.find("Material") == root.end())
    {
        throw std::domain_error("Material is not in input file");
    }
    if (root.find("SimulationCases") == root.end())
    {
        throw std::domain_error("SimulationCases is not in input file");
    }
}

std::unordered_set<std::string> simulation_parser::parse_material_names(json const& materials) const
{
    std::unordered_set<std::string> material_names;

    for (auto const& material : root["Material"])
    {
        if (material.find("Name") == material.end())
        {
            throw std::domain_error("Material: Name is missing");
        }

        auto const [it, inserted] = material_names.emplace(material["Name"].get<std::string>());

        if (!inserted) throw std::domain_error("Material was not parsed correctly");
    }
    return material_names;
}

std::unordered_set<std::string> simulation_parser::parse_part_names(
    json const& parts,
    std::unordered_set<std::string> const& material_names) const
{
    std::unordered_set<std::string> part_names;

    // Load in all the part names and error check
    for (auto const& part : root["Part"])
    {
        if (part.find("Name") == part.end())
        {
            throw std::domain_error("Part: Name is missing");
        }

        if (ranges::find(material_names, part["Material"].get<std::string>()) == material_names.end())
        {
            throw std::domain_error("The part material was not found in the provided "
                                    "materials\n");
        }

        auto const [it, inserted] = part_names.emplace(part["Name"].get<std::string>());

        if (!inserted) throw std::domain_error("Part is defined more than once");
    }
    return part_names;
}
}
