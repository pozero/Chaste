#include "SimplePolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"

SimplePolarityModifier::SimplePolarityModifier(double polarity_delta, std::string polarity_name) : mPolarityDelta(polarity_delta), mPolarityName(std::move(polarity_name))
{
}

SimplePolarityModifier::~SimplePolarityModifier()
{
}

void SimplePolarityModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double theta = p_cell->GetCellData()->GetItem(mPolarityName);
        double const delta_theta = 2.0 * (p_gen->ranf() - 0.5) * mPolarityDelta;
        theta += delta_theta;
        p_cell->GetCellData()->SetItem(mPolarityName, theta);
    }
}

void SimplePolarityModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double const theta = 2.0 * M_PI * p_gen->ranf();
        p_cell->GetCellData()->SetItem(mPolarityName, theta);
    }
}

void SimplePolarityModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}
