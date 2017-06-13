/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 * Copyright Darcy Beurle, 2016.
 */

#pragma once

#include "PlasticMaterial.hpp"

namespace neon
{
class PerfectPlasticElastic : public PlasticMaterial
{
public:
    PerfectPlasticElastic(Json::Value const& material_data);

    virtual ~PerfectPlasticElastic() = default;

    /**
     * Return the yield stress without to no hardening
     * @param effective_strain - The Von Mises effective strain
     * @return the stress
     */
    virtual double evaluate_yield_function(double effective_strain) const override;

    /** @return plastic modulus of a EPP material (0.0) */
    virtual double plastic_modulus(double effective_strain) const override;

protected:
    double yield_stress;
};
}
