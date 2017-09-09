
#pragma once

#include "AbstractModule.hpp"

#include "assembler/solid/femStaticMatrix.hpp"
#include "mesh/solid/femMesh.hpp"

namespace neon
{
/**
 * SolidMechanicsModule is responsible for handling the setup and simulation of the class
 * of three dimensional solid mechanics problems.
 */
class SolidMechanicsModule : public AbstractModule
{
public:
    SolidMechanicsModule(BasicMesh const& mesh,
                         Json::Value const& material,
                         Json::Value const& simulation);

    virtual void perform_simulation() override final { fem_matrix.solve(); }

protected:
    solid::femMesh fem_mesh;           //!< Mesh with the solid routines
    solid::femStaticMatrix fem_matrix; //!< Nonlinear solver routines
};
}
