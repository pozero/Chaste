#include "BoxBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <algorithm>

double clamp(double v, double min, double max)
{
    if (v < min)
    {
        return min;
    }
    else if (v > max)
    {
        return max;
    }
    return v;
}

BoxBoundaryCondition::BoxBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation, double top, double bottom, double left, double right)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation), m_top(top), m_bottom(bottom), m_left(left), m_right(right)
{
    if (dynamic_cast<VertexBasedCellPopulation<2>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("A VertexBasedCellPopulation must be used with the box boundary condition object");
    }
}

void BoxBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    unsigned const num_node = mpCellPopulation->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_node; ++node_idx)
    {
        Node<2>* node = mpCellPopulation->GetNode(node_idx);
        c_vector<double, 2>& node_location = node->rGetModifiableLocation();
        node_location[0] = clamp(node_location[0], m_left, m_right);
        node_location[1] = clamp(node_location[1], m_bottom, m_top);
    }
}

bool BoxBoundaryCondition::VerifyBoundaryCondition()
{
    bool satisfied = true;
    unsigned const num_node = mpCellPopulation->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_node; ++node_idx)
    {
        Node<2>* node = mpCellPopulation->GetNode(node_idx);
        c_vector<double, 2> const& node_location = node->rGetLocation();
        if (node_location[0] < m_left || node_location[0] > m_right || node_location[1] < m_bottom || node_location[1] > m_top)
        {
            satisfied = false;
            break;
        }
    }
    return satisfied;
}

void BoxBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
}
