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

public:
    PolarityModifier();

    virtual ~PolarityModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double GetNeighborAlignmentIntensity() const;

    double GetShapeAlignmentIntensity() const;

    double GetProtrusionAlignmentIntensity() const;

    void SetNeighborAlignmentIntensity(double newParam);

    void SetShapeAlignmentIntensity(double newParam);

    void SetProtrusionAlignmentIntensity(double newParam);

    void SetDt(double dt);
};

#endif
