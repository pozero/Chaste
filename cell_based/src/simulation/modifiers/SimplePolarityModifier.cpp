#include "SimplePolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"

template <unsigned DIM>
SimplePolarityModifier<DIM>::SimplePolarityModifier()
        : mPolarityDelta(1.0)
{
}

template <unsigned DIM>
SimplePolarityModifier<DIM>::~SimplePolarityModifier()
{
}

template <unsigned DIM>
void SimplePolarityModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
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

template <unsigned DIM>
void SimplePolarityModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double const theta = 2.0 * M_PI * p_gen->ranf();
        p_cell->GetCellData()->SetItem("polarity theta", theta);
    }
}

template <unsigned DIM>
void SimplePolarityModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PolarityDelta>" << mPolarityDelta << "</PolarityDelta>\n";
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template <unsigned DIM>
double SimplePolarityModifier<DIM>::GetPolarityDelta() const
{
    return mPolarityDelta;
}

template <unsigned DIM>
void SimplePolarityModifier<DIM>::SetPolarityDelta(double newParam)
{
    mPolarityDelta = newParam;
}

template class SimplePolarityModifier<1>;
template class SimplePolarityModifier<2>;
template class SimplePolarityModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimplePolarityModifier)
