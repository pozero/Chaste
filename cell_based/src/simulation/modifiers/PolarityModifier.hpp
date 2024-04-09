#ifndef POLARITYMODIFIER_HPP_
#define POLARITYMODIFIER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"

class PolarityModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mNeighborAlignmentIntensity;

    double mShapeAlignmentIntensity;

    double mProtrusionAlignmentIntensity;

    double mDt;

    std::string mPolarityName;

public:
    PolarityModifier(double neighbor, double shape, double protrusion, double dt, std::string polarity_name);

    virtual ~PolarityModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#endif
