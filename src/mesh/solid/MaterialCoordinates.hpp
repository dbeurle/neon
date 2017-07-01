
#pragma once

#include "mesh/NodalCoordinates.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace neon::solid
{
class MaterialCoordinates : public NodalCoordinates
{
public:
    /** Construct this class using a set of initial coordinates */
    MaterialCoordinates(Vector const& initial_coordinates);

    /** @return element reference configuration based on the local node numbers*/
    Matrix initial_configuration(List const& local_nodes) const;

    /** @return element current configuration based on the local node numbers*/
    Matrix current_configuration(List const& local_nodes) const;

    /** @param du incremental displacement vector (x,y,z...) */
    void update_current_configuration(Vector const& du) { x += du; };

    Vector displacement() const { return x - X; }

    Vector displacement(List const& local_dofs) const;

protected:
    Matrix get_configuration(List const& local_nodes, Vector const& configuration) const;

protected:
    Vector x; //!< Current configuration
};

inline Matrix MaterialCoordinates::initial_configuration(List const& local_node_list) const
{
    return this->get_configuration(local_node_list, X);
}

inline Matrix MaterialCoordinates::current_configuration(List const& local_node_list) const
{
    return this->get_configuration(local_node_list, x);
}
}
