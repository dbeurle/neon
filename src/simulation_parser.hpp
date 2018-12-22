
#pragma once

/// @file

#include <list>
#include <map>
#include <memory>
#include <string>
#include <set>
#include <utility>

#include "mesh/basic_mesh.hpp"
#include "io/json_forward.hpp"

namespace neon
{
// Forward declarations
namespace geometry
{
class profile;
}
class abstract_module;

class simulation_parser
{
public:
    explicit simulation_parser(std::string const& input_file_name);

    ~simulation_parser();

    void start();

    /// Number of hardware threads to use
    static int threads;

protected:
    void parse();

    void build_simulation_tree();

    /// Find all of the children that depend on a given simulation case
    void find_children(std::string const& parent_name, std::string const& next_parent_name);

    /// Print the welcome banner to the terminal
    void print_banner() const;

    void check_input_fields() const;

    /// Extract the material names from the input file
    [[nodiscard]] std::set<std::string> allocate_material_names(json const& materials) const;

    /// Extract the part names from the input file
    [[nodiscard]] std::set<std::string> parse_part_names(
        json const& parts,
        std::set<std::string> const& material_names) const;

    [[nodiscard]] std::set<std::string> allocate_profile_names(json const& profiles);

protected:
    /// Input data file name
    std::string input_file_name;

    /// Store the name as key and the mesh and material as a value pair
    std::map<std::string, std::pair<basic_mesh, json>> mesh_store;
    /// Profiles from input file root
    std::map<std::string, std::unique_ptr<geometry::profile>> profile_store;

    std::map<std::string, std::list<json>> multistep_simulations;

    /// Simulation modules to run
    std::vector<std::unique_ptr<abstract_module>> modules;

    /// The file input
    json root;
};
}
