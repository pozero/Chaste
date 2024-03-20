#ifndef FIXEDTARGETAREAMODIFIER_HPP_
#define FIXEDTARGETAREAMODIFIER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "TruncatedNormal.hpp"

#include <vector>

#define FIXED_TARGET_AREA_SIZE_TAG "fixed target area size tag"

#define FIXED_TARGET_AREA_SMALL 1.0
#define FIXED_TARGET_AREA_MEDIUM 2.0
#define FIXED_TARGET_AREA_LARGE 3.0

class FixedTargetAreaModifier : public AbstractCellBasedSimulationModifier<2, 2>
{
private:
    TruncatedNormal const* mTruncatedNormal = nullptr;

    double mShapeIndex;

    double mAverageTargetArea;

    double mAverageTargetPerimeter;

public:
    FixedTargetAreaModifier(TruncatedNormal const& truncated_normal, double shape_index);

    FixedTargetAreaModifier(double uniform_target_area, double shape_index);

    virtual ~FixedTargetAreaModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2, 2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2, 2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#endif /*FIXEDTARGETAREAMODIFIER_HPP_*/
