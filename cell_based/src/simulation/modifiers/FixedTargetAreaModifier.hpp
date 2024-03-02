#ifndef FIXEDTARGETAREAMODIFIER_HPP_
#define FIXEDTARGETAREAMODIFIER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"

#include <vector>

#define FIXED_TARGET_AREA_SIZE_TAG "fixed target area size tag"

#define FIXED_TARGET_AREA_SMALL 1.0
#define FIXED_TARGET_AREA_MEDIUM 2.0
#define FIXED_TARGET_AREA_LARGE 3.0

template <unsigned DIM>
class FixedTargetAreaModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:
    std::vector<double> mCellTargetArea;

    std::vector<double> mCellTargetPerimeter;

    std::vector<double> mCellTargetAreaCategory;

    double mParentNormalMean;

    double mParentNormalVariance;

    double mShapeIndex;

    double TruncatedNormal() const;

    double ParentNormalCDF(double val) const;

    double TruncatedNormalQuantile(double phi_a, double phi_b, double precentage) const;

public:
    FixedTargetAreaModifier(unsigned cell_count, double parent_normal_mean, double parent_normal_variance, double shape_index);

    FixedTargetAreaModifier(unsigned cell_count, double uniform_target_area, double shape_index);

    virtual ~FixedTargetAreaModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#endif /*FIXEDTARGETAREAMODIFIER_HPP_*/
