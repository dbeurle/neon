
#pragma once

#include "mesh/BasicMesh.hpp"

#include <map>
#include <string>
#include <utility>

#include <json/value.h>

namespace neon
{
class SimulationControl
{
public:
    explicit SimulationControl(std::string const& input_file_name);

    auto number_of_threads() const { return threads; }

    void start();

protected:
    void parse();

    void build_simulation_tree();

    /** Find all of the children that depend on a given simulation case */
    void find_children(std::string const& parent_name, std::string const& next_parent_name);

protected:
    int threads = 1; //!< Number of hardware threads to use

    std::string input_file_name;

    // Store the name, mesh connectivity and material
    std::map<std::string, std::pair<BasicMesh, Json::Value>> mesh_store;

    std::map<std::string, std::vector<Json::Value>> multistep_simulations;

    Json::Value root; // The file input
};
}
