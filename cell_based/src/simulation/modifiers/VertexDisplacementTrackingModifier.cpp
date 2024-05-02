#include "VertexDisplacementTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

VertexDisplacementTrackingModifier::VertexDisplacementTrackingModifier(double dt, std::string const& output_file)
        : AbstractCellBasedSimulationModifier<2>(), mDt(dt), mOutputFileName(output_file), mMSD(1, 0.0)
{
}

VertexDisplacementTrackingModifier::~VertexDisplacementTrackingModifier()
{
    mOutputFile << "Time, Displacement\n";
    for (uint32_t i = 0; i < mMSD.size(); ++i)
    {
        double const time = i * mDt;
        double const displacement = mMSD[i];
        mOutputFile << time << ", " << displacement << '\n';
    }
    mOutputFile.close();
}

void VertexDisplacementTrackingModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    VertexBasedCellPopulation<2>* cell_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_cells = cell_population->GetNumElements();
    double msd = 0.0;
    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        auto cell_data = cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData();
        double const centroid_x = cell_data->GetItem("centroid x");
        double const centroid_y = cell_data->GetItem("centroid y");
        double const initial_centroid_x = cell_data->GetItem("initial centroid x");
        double const initial_centroid_y = cell_data->GetItem("initial centroid y");
        double const squared_displacement = (centroid_x - initial_centroid_x) * (centroid_x - initial_centroid_x) + (centroid_y - initial_centroid_y) * (centroid_y - initial_centroid_y);
        msd += squared_displacement;
    }
    msd /= num_cells;
    if (std::isnan(msd))
    {
        EXCEPTION("MSD isn't a number");
    }
    mMSD.push_back(msd);
}

void VertexDisplacementTrackingModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory)
{
    FileFinder file_finder{ mOutputFileName + ".csv", RelativeTo::ChasteTestOutput };
    mOutputFile.open(file_finder.GetAbsolutePath());
    if (!mOutputFile.is_open())
    {
        EXCEPTION("Can't open file " << file_finder.GetAbsolutePath());
    }
}

void VertexDisplacementTrackingModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}
