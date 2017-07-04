
#pragma once

#include "mesh/BasicMesh.hpp"

#include "assembler/solid/femMatrix.hpp"
#include "mesh/solid/femMesh.hpp"

#include <map>
#include <string>
#include <utility>

namespace neon
{
class SimulationControl
{
public:
    SimulationControl(std::string const& input_file_name);

    auto number_of_threads() const { return threads; }

    void start();

protected:
    void parse();

protected:
    int threads = 1; //!< Number of hardware threads to use

    std::string input_file_name;

    // Store the name, material and mesh connectivity
    std::map<std::string, BasicMesh> mesh_store;
};
}
