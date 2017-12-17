
#pragma once

#include "mesh/NodalCoordinates.hpp"

#include "vtkSmartPointer.h"
class vtkPoints;
class vtkDoubleArray;

namespace neon
{
class MaterialCoordinates : public NodalCoordinates
{
public:
    /** Construct this class using a set of initial coordinates */
    MaterialCoordinates(vector const& initial_coordinates);

    /** @return element reference configuration based on the local node numbers*/
    [[nodiscard]] Matrix3x initial_configuration(local_indices const& local_nodes) const;

    /** @return element current configuration based on the local node numbers*/
    [[nodiscard]] Matrix3x current_configuration(local_indices const& local_nodes) const;

    /** @param u - displacement vector from initial configuration (x,y,z...) */
    void update_current_configuration(vector const& u) { x = X + u; };

    /** @param u - displacement vector from initial configuration (x,y...) */
    void update_current_xy_configuration(vector const& u);

    [[nodiscard]] Vector displacement() const { return x - X; }

    [[nodiscard]] Vector displacement(local_indices const& local_dofs) const;

    /** @return a vtk object of the initial coordinates */
    [[nodiscard]] vtkSmartPointer<vtkPoints> vtk_coordinates() const;

    /** @return a vtk array of nodal displacements */
    [[nodiscard]] vtkSmartPointer<vtkDoubleArray> vtk_displacement() const;

protected:
    [[nodiscard]] Matrix3x get_configuration(local_indices const& local_nodes,
                                             vector const& configuration) const;

protected:
    Vector x; //!< Current configuration
};

inline Matrix3x MaterialCoordinates::initial_configuration(List const& local_node_list) const
{
    return this->get_configuration(local_node_list, X);
}

inline Matrix3x MaterialCoordinates::current_configuration(List const& local_node_list) const
{
    return this->get_configuration(local_node_list, x);
}
}
