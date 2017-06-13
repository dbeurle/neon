
#pragma once

#include "PerfectPlasticElastic.hpp"

namespace neon
{
class IsotropicPlasticElastic : public PerfectPlasticElastic
{
public:
    IsotropicPlasticElastic(const Json::Value& materialData);

    virtual ~IsotropicPlasticElastic() = default;

    /**
     * Return the yield stress due to no hardening
     * @param effective_strain
     * @return
     */
    virtual double evaluate_yield_function(double effective_strain) const override;

    /**
     * Return the plastic modulus of a EPP material which is zero
     * @param effective_strain
     * @return
     */
    virtual double plastic_modulus(double effective_strain) const override;

protected:
    double failure_stress; //!< Material ultimate tensile stress
    double failure_strain; //!< Material ultimate tensile strain
};
}
