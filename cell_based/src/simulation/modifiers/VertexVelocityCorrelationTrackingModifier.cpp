#include "VertexVelocityCorrelationTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

double get_distance(c_vector<double, 2> const& v1, c_vector<double, 2> const& v2)
{
    c_vector<double, 2> const diff = v1 - v2;
    return std::sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
}

double get_correlation_cofficient(c_vector<double, 2> const& v1, c_vector<double, 2> const& v2)
{
    double const dot_product = v1[0] * v2[0] + v1[1] * v2[1];
    double const norm1 = v1[0] * v1[0] + v1[1] * v1[1];
    double const norm2 = v2[0] * v2[0] + v2[1] * v2[1];
    double const norm_product = std::sqrt(norm1 * norm2);
    return dot_product / norm_product;
}

VertexVelocityCorrelationTrackingModifier::VertexVelocityCorrelationTrackingModifier(double bucket_width, std::string const& output_file)
        : AbstractCellBasedSimulationModifier<2>(), mBucketWidth(bucket_width), mOutputFileName(output_file)
{
}

VertexVelocityCorrelationTrackingModifier::~VertexVelocityCorrelationTrackingModifier()
{
    mOutputFile << "Distance, Correlation\n";
    for (auto iter = mCorrelationCoefficients.begin(); iter != mCorrelationCoefficients.end(); ++iter)
    {
        double const distance = iter->first * mBucketWidth;
        int const count = mCorrelationCoefficientCounts.find(iter->first)->second;
        double const average_correlation = iter->second / static_cast<double>(count);
        mOutputFile << distance << ", " << average_correlation << '\n';
    }
    mOutputFile.close();
}

void VertexVelocityCorrelationTrackingModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    VertexBasedCellPopulation<2>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_nodes = p_cell_population->GetNumNodes();
    std::vector<c_vector<double, 2> > node_velocities(num_nodes);
    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        Node<2>* node = p_cell_population->GetNode(node_idx);
        std::vector<double>& node_attributes = node->rGetNodeAttributes();
        node_velocities[node_idx][0] = node_attributes[2];
        node_velocities[node_idx][1] = node_attributes[3];
    }

    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const node_cnt = elem_iter->GetNumNodes();
        for (unsigned i = 0; i < node_cnt; ++i)
        {
            unsigned const next_i = (i + 1) % node_cnt;
            unsigned const node_idx = elem_iter->GetNode(i)->GetIndex();
            unsigned const next_node_idx = elem_iter->GetNode(next_i)->GetIndex();
            c_vector<double, 2> const node_velocity = node_velocities[node_idx];
            c_vector<double, 2> const next_node_velocity = node_velocities[next_node_idx];
            double const veloctiy_correlation_cofficient = get_correlation_cofficient(node_velocity, next_node_velocity);
            c_vector<double, 2> const node_location = elem_iter->GetNodeLocation(i);
            c_vector<double, 2> const next_node_location = elem_iter->GetNodeLocation(next_i);
            double const distance = get_distance(node_location, next_node_location);
            int const bucket_idx = std::floor(distance / mBucketWidth);
            mCorrelationCoefficients[bucket_idx] += veloctiy_correlation_cofficient;
            mCorrelationCoefficientCounts[bucket_idx] += 1;
        }
    }
}

void VertexVelocityCorrelationTrackingModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory)
{
    FileFinder file_finder{ mOutputFileName + ".csv", RelativeTo::ChasteTestOutput };
    mOutputFile.open(file_finder.GetAbsolutePath());
    if (!mOutputFile.is_open())
    {
        EXCEPTION("Can't open file " << file_finder.GetAbsolutePath());
    }
}

void VertexVelocityCorrelationTrackingModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}
