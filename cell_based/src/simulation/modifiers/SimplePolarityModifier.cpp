#include "SimplePolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"

SimplePolarityModifier::SimplePolarityModifier()
        : mPolarityDelta(1.0)
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
        double theta = p_cell->GetCellData()->GetItem("polarity theta");
        double const delta_theta = 2.0 * (p_gen->ranf() - 0.5) * mPolarityDelta;
        theta += delta_theta;
        p_cell->GetCellData()->SetItem("polarity theta", theta);
    }
}

void SimplePolarityModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double const theta = 2.0 * M_PI * p_gen->ranf();
        p_cell->GetCellData()->SetItem("polarity theta", theta);
    }
}

void SimplePolarityModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PolarityDelta>" << mPolarityDelta << "</PolarityDelta>\n";
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

double SimplePolarityModifier::GetPolarityDelta() const
{
    return mPolarityDelta;
}

void SimplePolarityModifier::SetPolarityDelta(double newParam)
{
    mPolarityDelta = newParam;
}
