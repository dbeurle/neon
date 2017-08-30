
#pragma once

#include "mesh/BasicMesh.hpp"
// #include "modules/AbstractModule.hpp"

#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include <json/value.h>

namespace neon
{
class SimulationControl
{
public:
    explicit SimulationControl(std::string const& input_file_name);

    void start(bool const is_initial_pass = true);

    static int threads; //!< Number of hardware threads to use
protected:
    void parse();

    void build_simulation_tree();

    /** Find all of the children that depend on a given simulation case */
    void find_children(std::string const& parent_name, std::string const& next_parent_name);

    /** Print the welcome banner to the terminal */
    void print_banner() const;

    /** Factory method to create simulation modules */
    std::unique_ptr<AbstractModule> make_module() const;

    void check_input_fields(Json::Value const& root) const;

    /** Extract the material names from the input file */
    std::unordered_set<std::string> parse_material_names(Json::Value const& materials) const;

    /** Extract the part names from the input file */
    std::unordered_set<std::string> parse_part_names(
        Json::Value const& parts, std::unordered_set<std::string> const& material_names) const;

protected:
    std::string input_file_name;

    // Store the name, mesh connectivity and material
    std::map<std::string, std::pair<BasicMesh, Json::Value>> mesh_store;

    std::map<std::string, std::vector<Json::Value>> multistep_simulations;

    // std::vector<std::unique_ptr<AbstractModule>> modules;

    Json::Value root; // The file input
};
}
