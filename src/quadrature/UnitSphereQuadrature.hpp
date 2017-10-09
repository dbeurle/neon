
#include "NumericalQuadrature.hpp"

namespace neon
{
/** UnitSphereQuadrature implements the symmetric quadrature rules */
class UnitSphereQuadrature : public VolumeQuadrature
{
public:
    /**
     * Quadrature nomenclature from
     * Numerical integration on the sphere and its effect on the material
     * symmetry of constitutive equations - A comparative study
     * by Ehret et.al.
     */
    enum class Rule { BO21, BO33 };

public:
    /** Fill the quadrature coordinates and weightings */
    UnitSphereQuadrature(Rule const rule = Rule::BO21);

protected:
    void precompute_coordinates();
};
}
