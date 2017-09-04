
#pragma once

#include "AbstractModule.hpp"

#include "mesh/solid/femMesh.hpp"

namespace neon
{
/**
 * SolidModule is responsible for handling the setup and simulation of the class
 * of three dimensional solid mechanics problems.
 */
class SolidModule : public AbstractModule
{
public:
    template <typename... args>
    SolidModule(args...) : fem_mesh(args...), fem_matrix()
    {
    }

    virtual void perform_simulation() override final {}

protected:
    solid::femMesh fem_mesh;     //!< Mesh with the solid routines
    solid::femMatrix fem_matrix; //!< Nonlinear solver routines
};
}
