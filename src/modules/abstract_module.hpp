
#pragma once

/// @file

#include "io/json.hpp"

namespace neon
{
/**
 * abstract_module is an abstract base class for the module system.  Derived modules
 * must provide methods that perform the solution routines, in addition to
 * apriori checks if the simulation is likely to complete and no errors are
 * found in the input file.
 */
class abstract_module
{
public:
    virtual void perform_simulation() = 0;

    virtual ~abstract_module() = default;

protected:
    std::map<std::string, std::vector<json>> multistep_simulations;
};
}
