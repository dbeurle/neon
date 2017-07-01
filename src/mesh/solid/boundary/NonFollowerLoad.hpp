
#pragma once

namespace neon::solid
{
/**
 * NonfollowerLoad is a base class for the class of loads which do not follow
 * the body throughout the simulation.  Every nonfollower load will contribute
 * to the external force vector, whether from a volume or a surface load
 */
class NonFollowerLoad
{
public:
    /** @return an element external force vector for a given element */
    virtual Vector external_force(int element) const = 0;

    auto elements() const { return nodal_connectivity.size(); }

protected:
    std::vector<List> nodal_connectivity;
    std::vector<List> dof_list;

    double prescribed_value;
};

class Traction : public NonFollowerLoad
{
public:
    Traction(std::unique_ptr<SurfaceInterpolation>&& sf) : sf(sf) {}

    virtual Vector external_force(int element) const override;

protected:
    std::unique_ptr<SurfaceInterpolation> sf;

private:
};
}
