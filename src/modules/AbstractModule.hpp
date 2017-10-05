
#pragma once

#include <json/value.h>

namespace neon
{
/**
 * AbstractModule is an abstract base class for the module system.  Derived modules
 * must provide methods that perform the solution routines, in addition to
 * apriori checks if the simulation is likely to complete and no errors are
 * found in the input file.
 */
class AbstractModule
{
public:
    virtual void perform_simulation() = 0;

    virtual ~AbstractModule() = default;

protected:
    std::map<std::string, std::vector<Json::Value>> multistep_simulations;
};
}
