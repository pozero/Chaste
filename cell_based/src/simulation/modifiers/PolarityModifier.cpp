#include "PolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"

PolarityModifier::PolarityModifier() : mNeighborAlignmentIntensity(0.05),
                                       mShapeAlignmentIntensity(0.1),
                                       // mProtrusionAlignmentIntensity(3.5),
                                       mDt(0.0)
{
}

PolarityModifier::~PolarityModifier()
{
}

void PolarityModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    double const one_dt = 1.0 / mDt;
    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    std::vector<double> velocity_thetas(p_cell_population->GetNumElements());
    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        boost::shared_ptr<CellData> p_cell_data = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData();
        unsigned const node_cnt = elem_iter->GetNumNodes();
        c_vector<double, 2> centroid = zero_vector<double>(2);
        for (unsigned i = 0; i < node_cnt; ++i)
        {
            centroid += elem_iter->GetNodeLocation(i);
        }
        centroid *= 1.0 / static_cast<double>(node_cnt);
        c_vector<double, 2> last_centroid = zero_vector<double>(2);
        last_centroid[0] = p_cell_data->GetItem("centroid x");
        last_centroid[1] = p_cell_data->GetItem("centroid y");
        c_vector<double, 2> const velocity = (centroid - last_centroid) * one_dt;
        p_cell_data->SetItem("centroid x", centroid[0]);
        p_cell_data->SetItem("centroid y", centroid[1]);
        p_cell_data->SetItem("velocity x", velocity[0]);
        p_cell_data->SetItem("velocity y", velocity[1]);
        velocity_thetas[elem_idx] = acos(velocity[0] / sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]));
    }
    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        double neighbor_alignment = 0.0;
        double shape_alignment = 0.0;
        // double protrusion_alignment = 0.0;
        unsigned const elem_idx = elem_iter->GetIndex();
        c_matrix<double, 2, 2> shape_tensor = zero_matrix<double>(2, 2);
        boost::shared_ptr<CellData> p_cell_data = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData();
        double theta = p_cell_data->GetItem("polarity theta");
        std::set<unsigned> neighbor_elem_indices{};
        unsigned const node_cnt = elem_iter->GetNumNodes();
        c_vector<double, 2> node_location = elem_iter->GetNodeLocation(0);
        for (unsigned i = 0; i < node_cnt; ++i)
        {
            std::set<unsigned> const& relavant_node_idx = elem_iter->GetNode(i)->rGetContainingElementIndices();
            neighbor_elem_indices.insert(relavant_node_idx.begin(), relavant_node_idx.end());

            unsigned next_i = (i + 1) % node_cnt;
            c_vector<double, 2> const next_node_location = elem_iter->GetNodeLocation(next_i);
            shape_tensor(0, 0) += (node_location[0] - next_node_location[0]) * (node_location[0] - next_node_location[0]);
            shape_tensor(0, 1) += (node_location[0] - next_node_location[0]) * (node_location[1] - next_node_location[1]);
            shape_tensor(1, 0) += (node_location[0] - next_node_location[0]) * (node_location[1] - next_node_location[1]);
            shape_tensor(1, 1) += (node_location[1] - next_node_location[1]) * (node_location[1] - next_node_location[1]);
            node_location = next_node_location;
        }

        neighbor_elem_indices.erase(elem_idx);
        for (unsigned const neighbor_elem_idx : neighbor_elem_indices)
        {
            neighbor_alignment += sin(velocity_thetas[neighbor_elem_idx] - theta);
        }

        shape_tensor *= 1.0 / static_cast<double>(node_cnt);
        double const eigenvector_theta_biggest_eigenvalue = [&shape_tensor]
        {
            double const diagonal_sum = shape_tensor(0, 0) + shape_tensor(1, 1);
            double const determinant = shape_tensor(0, 0) * shape_tensor(1, 1) - shape_tensor(0, 1) * shape_tensor(1, 0);
            double const lambda_biggest = 0.5 * diagonal_sum + sqrt(0.25 * diagonal_sum * diagonal_sum - determinant);
            if (abs(shape_tensor(1, 0)) >= 1e-10)
            {
                return acos((lambda_biggest - shape_tensor(1, 1)) / sqrt((lambda_biggest - shape_tensor(1, 1))) * (lambda_biggest - shape_tensor(1, 1)) + shape_tensor(1, 0) * shape_tensor(1, 0));
            }
            else if (abs(shape_tensor(0, 1)) >= 1e-10)
            {
                return acos(shape_tensor(0, 1) / sqrt(shape_tensor(0, 1) * shape_tensor(0, 1) + (lambda_biggest - shape_tensor(0, 0)) * (lambda_biggest - shape_tensor(0, 0))));
            }
            return 0.0;
        }();
        shape_alignment = sin(eigenvector_theta_biggest_eigenvalue - theta);

        // theta += mDt * (mNeighborAlignmentIntensity * neighbor_alignment + mShapeAlignmentIntensity * shape_alignment + mProtrusionAlignmentIntensity * protrusion_alignment);
        theta += mDt * (mNeighborAlignmentIntensity * neighbor_alignment + mShapeAlignmentIntensity * shape_alignment);
        p_cell_data->SetItem("polarity theta", theta);
    }
}

void PolarityModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory)
{
    if (mDt == 0.0)
    {
        EXCEPTION("Delta t must be set before simulation");
    }

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    for (CellPtr p_cell : rCellPopulation.rGetCells())
    {
        double const theta = 2.0 * M_PI * p_gen->ranf();
        p_cell->GetCellData()->SetItem("polarity theta", theta);
    }

    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    // mCellOutwradVectorSum.resize(p_cell_population->GetNumElements(), zero_vector<double>(2));
    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        unsigned const node_cnt = elem_iter->GetNumNodes();
        c_vector<double, 2> centroid = zero_vector<double>(2);
        for (unsigned i = 0; i < node_cnt; ++i)
        {
            centroid += elem_iter->GetNodeLocation(i);
        }
        centroid *= 1.0 / static_cast<double>(node_cnt);
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("centroid x", centroid[0]);
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("centroid y", centroid[1]);
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("velocity x", 0.0);
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("velocity y", 0.0);
    }
}

void PolarityModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NeighborAlignmentIntensity>" << mNeighborAlignmentIntensity << "</NeighborAlignmentIntensity>\n";
    *rParamsFile << "\t\t\t<ShapeAlignmentIntensity>" << mShapeAlignmentIntensity << "</ShapeAlignmentIntensity>\n";
    // *rParamsFile << "\t\t\t<ProtrusionAlignmentIntensity>" << mProtrusionAlignmentIntensity << "</ProtrusionAlignmentIntensity>\n";
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

double PolarityModifier::GetNeighborAlignmentIntensity() const
{
    return mNeighborAlignmentIntensity;
}
double PolarityModifier::GetShapeAlignmentIntensity() const
{
    return mShapeAlignmentIntensity;
}

// double PolarityModifier::GetProtrusionAlignmentIntensity() const
// {
// return mProtrusionAlignmentIntensity;
// }

void PolarityModifier::SetNeighborAlignmentIntensity(double newParam)
{
    mNeighborAlignmentIntensity = newParam;
}

void PolarityModifier::SetShapeAlignmentIntensity(double newParam)
{
    mShapeAlignmentIntensity = newParam;
}

// void PolarityModifier::SetProtrusionAlignmentIntensity(double newParam)
// {
// mProtrusionAlignmentIntensity = newParam;
// }

void PolarityModifier::SetDt(double dt)
{
    mDt = dt;
}
