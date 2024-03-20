#include "FixedTargetAreaModifier.hpp"
#include "RandomNumberGenerator.hpp"

#include <boost/math/special_functions/erf.hpp>

FixedTargetAreaModifier::FixedTargetAreaModifier(TruncatedNormal const& truncated_normal, double shape_index)
        : mTruncatedNormal(&truncated_normal),
          mShapeIndex(shape_index),
          mAverageTargetArea(truncated_normal.GetMean()),
          mAverageTargetPerimeter(mShapeIndex * sqrt(mAverageTargetArea))
{
}

FixedTargetAreaModifier::FixedTargetAreaModifier(double uniform_target_area, double shape_index)
        : mShapeIndex(shape_index),
          mAverageTargetArea(uniform_target_area),
          mAverageTargetPerimeter(mShapeIndex * sqrt(mAverageTargetArea))
{
}

FixedTargetAreaModifier::~FixedTargetAreaModifier()
{
}

void FixedTargetAreaModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2, 2>& rCellPopulation)
{
}

void FixedTargetAreaModifier::SetupSolve(AbstractCellPopulation<2, 2>& rCellPopulation, std::string outputDirectory)
{
    // nonuniform
    if (mTruncatedNormal)
    {
        double const phi_min = mTruncatedNormal->ParentNormalCDF(mTruncatedNormal->GetMin());
        double const phi_max = mTruncatedNormal->ParentNormalCDF(mTruncatedNormal->GetMax());
        double const quantile_33 = mTruncatedNormal->TruncatedNormalQuantile(phi_min, phi_max, 0.33);
        double const quantile_66 = mTruncatedNormal->TruncatedNormalQuantile(phi_min, phi_max, 0.66);
        std::vector<double> const& cell_target_areas = mTruncatedNormal->GetGenerated();
        std::vector<double> cell_target_perimeters{};
        std::vector<double> cell_area_tag{};
        cell_target_perimeters.reserve(cell_target_areas.size());
        cell_area_tag.reserve(cell_target_areas.size());
        for (double const target_area : cell_target_areas)
        {
            double const target_perimeter = mShapeIndex * std::sqrt(target_area);
            cell_target_perimeters.push_back(target_perimeter);
            if (target_area <= quantile_33)
            {
                cell_area_tag.push_back(FIXED_TARGET_AREA_SMALL);
            }
            else if (target_area <= quantile_66)
            {
                cell_area_tag.push_back(FIXED_TARGET_AREA_MEDIUM);
            }
            else /* target_area > quantile_66 */
            {
                cell_area_tag.push_back(FIXED_TARGET_AREA_LARGE);
            }
        }
        unsigned int idx = 0;
        for (CellPtr p_cell : rCellPopulation.rGetCells())
        {
            p_cell->GetCellData()->SetItem("target area", cell_target_areas[idx]);
            p_cell->GetCellData()->SetItem("target perimeter", cell_target_perimeters[idx]);
            p_cell->GetCellData()->SetItem(FIXED_TARGET_AREA_SIZE_TAG, cell_area_tag[idx]);
            ++idx;
        }
    }
    // uniform
    else
    {
        for (CellPtr p_cell : rCellPopulation.rGetCells())
        {
            p_cell->GetCellData()->SetItem("target area", mAverageTargetArea);
            p_cell->GetCellData()->SetItem("target perimeter", mAverageTargetPerimeter);
            p_cell->GetCellData()->SetItem(FIXED_TARGET_AREA_SIZE_TAG, FIXED_TARGET_AREA_MEDIUM);
        }
    }
}

void FixedTargetAreaModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}
