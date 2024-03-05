#include "PolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"

static double get_theta(c_vector<double, 2> const& vector)
{
    double const length = sqrt(vector[0] * vector[0] + vector[1] * vector[1]);
    c_vector<double, 2> const normalized = vector / length;
    return acos(normalized[0]);
}

static c_vector<double, 2> get_eigenvector_biggest_eigenvalue(c_matrix<double, 2, 2> const& shape_tensor)
{
    double const diagonal_sum = shape_tensor(0, 0) + shape_tensor(1, 1);
    double const determinant = shape_tensor(0, 0) * shape_tensor(1, 1) - shape_tensor(0, 1) * shape_tensor(1, 0);
    double const lambda_biggest = 0.5 * diagonal_sum + sqrt(0.25 * diagonal_sum * diagonal_sum - determinant);
    c_vector<double, 2> ret = zero_vector<double>(2);
    ret[0] = 1.0;
    ret[1] = 1.0;
    if (abs(shape_tensor(1, 0)) >= 1e-10)
    {
        ret[0] = lambda_biggest - shape_tensor(1, 1);
        ret[1] = shape_tensor(1, 0);
    }
    else if (abs(shape_tensor(0, 1)) >= 1e-10)
    {
        ret[0] = shape_tensor(0, 1);
        ret[1] = lambda_biggest - shape_tensor(0, 0);
    }
    return ret;
}

PolarityModifier::PolarityModifier() : mNeighborAlignmentIntensity(0.05),
                                       mShapeAlignmentIntensity(0.1),
                                       mProtrusionAlignmentIntensity(3.5),
                                       mDt(0.0)
{
}

PolarityModifier::~PolarityModifier()
{
}

void PolarityModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    double const one_dt = 1.0 / mDt;
    VertexBasedCellPopulation<2>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_cells = p_cell_population->GetNumElements();
    std::vector<double> velocity_thetas(num_cells);
    std::vector<c_vector<double, 2> > cell_last_centroid(num_cells, c_vector<double, 2>{});
    std::vector<c_vector<double, 2> > cell_centroid(num_cells, c_vector<double, 2>{});
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
        cell_centroid[elem_idx] = centroid;
        c_vector<double, 2> last_centroid = zero_vector<double>(2);
        last_centroid[0] = p_cell_data->GetItem("centroid x");
        last_centroid[1] = p_cell_data->GetItem("centroid y");
        cell_last_centroid[elem_idx] = last_centroid;
        c_vector<double, 2> const velocity = (centroid - last_centroid) * one_dt;
        p_cell_data->SetItem("centroid x", centroid[0]);
        p_cell_data->SetItem("centroid y", centroid[1]);
        p_cell_data->SetItem("velocity x", velocity[0]);
        p_cell_data->SetItem("velocity y", velocity[1]);
        velocity_thetas[elem_idx] = get_theta(velocity);
    }

    unsigned const num_nodes = p_cell_population->GetNumNodes();
    std::vector<c_vector<double, 2> > node_last_location(num_nodes, c_vector<double, 2>{});
    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        Node<2>* node = p_cell_population->GetNode(node_idx);
        c_vector<double, 2> const& current_location = node->rGetLocation();
        std::vector<double>& node_attributes = node->rGetNodeAttributes();
        if (node_attributes.size() == 4)
        {
            c_vector<double, 2> last_location = zero_vector<double>(2);
            last_location[0] = node_attributes[0];
            last_location[1] = node_attributes[1];
            node_last_location[node_idx] = last_location;
            c_vector<double, 2> const velocity = (current_location - last_location) * one_dt;
            node_attributes[0] = current_location[0];
            node_attributes[1] = current_location[1];
            node_attributes[2] = velocity[0];
            node_attributes[3] = velocity[1];
        }
        else if (node_attributes.empty())
        {
            // new node
            node_last_location[node_idx] = current_location;
            node->AddNodeAttribute(current_location[0]);
            node->AddNodeAttribute(current_location[1]);
            node->AddNodeAttribute(0.0);
            node->AddNodeAttribute(0.0);
        }
        else
        {
            EXCEPTION("Incomplete node with " << node_attributes.size() << " attributes found.");
        }
    }

    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        double neighbor_alignment = 0.0;
        double shape_alignment = 0.0;
        double protrusion_alignment = 0.0;
        unsigned const elem_idx = elem_iter->GetIndex();
        c_matrix<double, 2, 2> shape_tensor = zero_matrix<double>(2, 2);
        boost::shared_ptr<CellData> p_cell_data = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData();
        double theta = p_cell_data->GetItem("polarity theta");
        std::set<unsigned> neighbor_elem_indices{};
        unsigned const node_cnt = elem_iter->GetNumNodes();
        c_vector<double, 2> node_location = elem_iter->GetNodeLocation(0);
        double biggest_protrusion_theta = std::numeric_limits<double>::min();
        for (unsigned i = 0; i < node_cnt; ++i)
        {
            std::set<unsigned> const& relavant_elem_idx = elem_iter->GetNode(i)->rGetContainingElementIndices();
            neighbor_elem_indices.insert(relavant_elem_idx.begin(), relavant_elem_idx.end());

            unsigned next_i = (i + 1) % node_cnt;
            c_vector<double, 2> const next_node_location = elem_iter->GetNodeLocation(next_i);
            shape_tensor(0, 0) += (node_location[0] - next_node_location[0]) * (node_location[0] - next_node_location[0]);
            shape_tensor(0, 1) += (node_location[0] - next_node_location[0]) * (node_location[1] - next_node_location[1]);
            shape_tensor(1, 0) += (node_location[0] - next_node_location[0]) * (node_location[1] - next_node_location[1]);
            shape_tensor(1, 1) += (node_location[1] - next_node_location[1]) * (node_location[1] - next_node_location[1]);

            unsigned const node_global_idx = elem_iter->GetNode(i)->GetIndex();
            c_vector<double, 2> const outward_vector = node_location - cell_centroid[elem_idx];
            c_vector<double, 2> const last_outward_vector = node_last_location[node_global_idx] - cell_last_centroid[elem_idx];
            c_vector<double, 2> const difference = outward_vector - last_outward_vector;
            double const theta = get_theta(difference);
            biggest_protrusion_theta = std::max(theta, biggest_protrusion_theta);

            node_location = next_node_location;
        }
        protrusion_alignment = sin(biggest_protrusion_theta - theta);

        neighbor_elem_indices.erase(elem_idx);
        for (unsigned const neighbor_elem_idx : neighbor_elem_indices)
        {
            neighbor_alignment += sin(velocity_thetas[neighbor_elem_idx] - theta);
        }

        shape_tensor *= 1.0 / static_cast<double>(node_cnt);
        c_vector<double, 2> const eigenvector_biggest_eigenvalue = get_eigenvector_biggest_eigenvalue(shape_tensor);
        double const eigenvector_theta_biggest_eigenvalue = get_theta(eigenvector_biggest_eigenvalue);
        shape_alignment = sin(eigenvector_theta_biggest_eigenvalue - theta);

        theta += mDt * (mNeighborAlignmentIntensity * neighbor_alignment + mShapeAlignmentIntensity * shape_alignment + mProtrusionAlignmentIntensity * protrusion_alignment);
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
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("initial centroid x", centroid[0]);
        p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->SetItem("initial centroid y", centroid[1]);
    }

    unsigned const num_nodes = p_cell_population->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        Node<2>* node = p_cell_population->GetNode(node_idx);
        c_vector<double, 2> const& node_location = node->rGetLocation();
        // location
        node->AddNodeAttribute(node_location[0]);
        node->AddNodeAttribute(node_location[1]);
        // velocity
        node->AddNodeAttribute(0.0);
        node->AddNodeAttribute(0.0);
    }
}

void PolarityModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
}

double PolarityModifier::GetNeighborAlignmentIntensity() const
{
    return mNeighborAlignmentIntensity;
}
double PolarityModifier::GetShapeAlignmentIntensity() const
{
    return mShapeAlignmentIntensity;
}

double PolarityModifier::GetProtrusionAlignmentIntensity() const
{
    return mProtrusionAlignmentIntensity;
}

void PolarityModifier::SetNeighborAlignmentIntensity(double newParam)
{
    mNeighborAlignmentIntensity = newParam;
}

void PolarityModifier::SetShapeAlignmentIntensity(double newParam)
{
    mShapeAlignmentIntensity = newParam;
}

void PolarityModifier::SetProtrusionAlignmentIntensity(double newParam)
{
    mProtrusionAlignmentIntensity = newParam;
}

void PolarityModifier::SetDt(double dt)
{
    mDt = dt;
}
