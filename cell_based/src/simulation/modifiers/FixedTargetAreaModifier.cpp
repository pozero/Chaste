#include "FixedTargetAreaModifier.hpp"
#include "RandomNumberGenerator.hpp"

#include <boost/math/special_functions/erf.hpp>

template <unsigned DIM>
FixedTargetAreaModifier<DIM>::FixedTargetAreaModifier() : AbstractCellBasedSimulationModifier<DIM>(),
                                                          mParentNormalMean(1.0),
                                                          mParentNormalVariance(1.0)
{
}

template <unsigned DIM>
FixedTargetAreaModifier<DIM>::~FixedTargetAreaModifier()
{
}

template <unsigned DIM>
void FixedTargetAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
}

template <unsigned DIM>
double FixedTargetAreaModifier<DIM>::TruncatedNormal() const
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
    while (val <= 0.0)
    {
        val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
    }
    return val;
}

template <unsigned DIM>
double FixedTargetAreaModifier<DIM>::ParentNormalCDF(double val) const
{
    return 0.5 * (1.0 + boost::math::erf((val - mParentNormalMean) / (sqrt(2.0) * mParentNormalVariance)));
}

template <unsigned DIM>
double FixedTargetAreaModifier<DIM>::TruncatedNormalQuantile(double phi_a, double phi_b, double percentage) const
{
    return sqrt(2.0) * mParentNormalVariance * boost::math::erf_inv(2.0 * percentage * (phi_b - phi_a) + 2.0 * phi_a - 1.0) + mParentNormalMean;
}

template <unsigned DIM>
void FixedTargetAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    double const phi_0 = ParentNormalCDF(0.0);
    double const phi_pos_inf = 1.0;
    double const quantile_33 = TruncatedNormalQuantile(phi_0, phi_pos_inf, 0.33);
    double const quantile_66 = TruncatedNormalQuantile(phi_0, phi_pos_inf, 0.66);
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double const SHAPE_INDEX = sqrt(24.0 / sqrt(3.0));
        double const target_area = TruncatedNormal();
        double const target_perimeter = SHAPE_INDEX * sqrt(target_area);
        p_cell->GetCellData()->SetItem("target area", target_area);
        p_cell->GetCellData()->SetItem("target perimeter", target_perimeter);
        if (target_area <= quantile_33)
        {
            p_cell->GetCellData()->SetItem(FIXED_TARGET_AREA_SIZE_TAG, FIXED_TARGET_AREA_SMALL);
        }
        else if (target_area <= quantile_66)
        {
            p_cell->GetCellData()->SetItem(FIXED_TARGET_AREA_SIZE_TAG, FIXED_TARGET_AREA_MEDIUM);
        }
        else /* target_area > quantile_66 */
        {
            p_cell->GetCellData()->SetItem(FIXED_TARGET_AREA_SIZE_TAG, FIXED_TARGET_AREA_LARGE);
        }
    }
}

template <unsigned DIM>
void FixedTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ParentNormalMean>" << mParentNormalMean << "</ParentNormalMean>\n";
    *rParamsFile << "\t\t\t<ParentNormalVariance>" << mParentNormalVariance << "</ParentNormalVariance>\n";
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template <unsigned DIM>
double FixedTargetAreaModifier<DIM>::GetParentNormalMean() const
{
    return mParentNormalMean;
}

template <unsigned DIM>
double FixedTargetAreaModifier<DIM>::GetParentNormalVariance() const
{
    return mParentNormalVariance;
}

template <unsigned DIM>
void FixedTargetAreaModifier<DIM>::SetParentNormalMean(double newParam)
{
    mParentNormalMean = newParam;
}

template <unsigned DIM>
void FixedTargetAreaModifier<DIM>::SetParentNormalVariance(double newParam)
{
    mParentNormalVariance = newParam;
}

template class FixedTargetAreaModifier<1>;
template class FixedTargetAreaModifier<2>;
template class FixedTargetAreaModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedTargetAreaModifier);
