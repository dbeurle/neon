
#pragma once

#include "AbstractModule.hpp"

#include "assembler/solid/femStaticMatrix.hpp"
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
    SolidModule(args...)
    {
    }

    virtual void perform_simulation() override final {}

protected:
    // solid::femMesh fem_mesh;           //!< Mesh with the solid routines
    // solid::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};
}
