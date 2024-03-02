#include "VertexDisplacementTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

VertexDisplacementTrackingModifier::VertexDisplacementTrackingModifier(double dt, std::string const& output_file)
        : AbstractCellBasedSimulationModifier<2>(), mDt(dt), mOutputFileName(output_file), mVerticesDisplacementAverage(1, 0.0)
{
}

VertexDisplacementTrackingModifier::~VertexDisplacementTrackingModifier()
{
    mOutputFile << "Time, Displacement\n";
    for (uint32_t i = 0; i < mVerticesDisplacementAverage.size(); ++i)
    {
        double const time = i * mDt;
        double const displacement = mVerticesDisplacementAverage[i];
        mOutputFile << time << ", " << displacement << '\n';
    }
    mOutputFile.close();
}

void VertexDisplacementTrackingModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    VertexBasedCellPopulation<2>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_nodes = p_cell_population->GetNumNodes();
    double displacement_squared_sum = 0.0;
    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        Node<2>* node = p_cell_population->GetNode(node_idx);
        std::vector<double>& node_attributes = node->rGetNodeAttributes();
        double const displacement = node_attributes[4] * node_attributes[4] + node_attributes[5] * node_attributes[5];
        displacement_squared_sum += displacement;
    }
    displacement_squared_sum /= static_cast<double>(num_nodes);
    mVerticesDisplacementAverage.push_back(mVerticesDisplacementAverage.back() + displacement_squared_sum);
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
