#include "VertexVelocityCorrelationTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

double get_correlation(c_vector<double, 2> const& v1, c_vector<double, 2> const& v2)
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
    for (auto& pair : mCorrelations)
    {
        double const distance_floor = pair.first;
        double const accumulated_correlation = pair.second;
        double const accumulated_distance = mDistances.find(distance_floor)->second;
        int const count = mBucketCounts.find(distance_floor)->second;
        double const distance = accumulated_distance / count;
        double const correlation = accumulated_correlation / count;
        mOutputFile << distance << ", " << correlation << "\n";
    }
    mOutputFile.close();
}

void VertexVelocityCorrelationTrackingModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    VertexBasedCellPopulation<2>* cell_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_cells = cell_population->GetNumAllCells();
    std::vector<c_vector<double, 2> > cell_velocities(num_cells);
    std::vector<c_vector<double, 2> > cell_centroids(num_cells);
    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        auto cell_data = cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData();
        cell_velocities[elem_idx][0] = cell_data->GetItem("velocity x");
        cell_velocities[elem_idx][1] = cell_data->GetItem("velocity y");
        cell_centroids[elem_idx][0] = cell_data->GetItem("centroid x");
        cell_centroids[elem_idx][1] = cell_data->GetItem("centroid y");
    }

    for (unsigned i = 0; i < num_cells; ++i)
    {
        auto const cell_i_velocity = cell_velocities[i];
        auto const cell_i_centroid = cell_centroids[i];
        for (unsigned j = 0; j < num_cells; ++j)
        {
            auto const cell_j_velocity = cell_velocities[j];
            auto const cell_j_centroid = cell_centroids[j];
            double const velocity_correlation = get_correlation(cell_i_velocity, cell_j_velocity);
            if (std::isnan(velocity_correlation))
            {
                continue;
            }
            double const distance = norm_2(cell_i_centroid - cell_j_centroid);
            double const distance_floor = std::floor(distance / mBucketWidth) * mBucketWidth;
            mCorrelations[distance_floor] += velocity_correlation;
            mDistances[distance_floor] += distance;
            mBucketCounts[distance_floor] += 1;
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
