#include "InelasticalCircleBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"

static double norm(c_vector<double, 2> const& vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
}

static c_vector<double, 2> normalize(c_vector<double, 2> const& vec)
{
    return vec / norm(vec);
}

InelasticalCircleBoundaryCondition::InelasticalCircleBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation,
                                                                       c_vector<double, 2> const& center,
                                                                       double radius) : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
                                                                                        m_circle_center(center),
                                                                                        m_circle_radius(radius)
{
    assert(radius > 0.0);
    if (dynamic_cast<VertexBasedCellPopulation<2>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("A VertexBasedCellPopulation must be used with the circle boundary condition object");
    }
}

void InelasticalCircleBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    unsigned const num_node = mpCellPopulation->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_node; ++node_idx)
    {
        Node<2>* node = mpCellPopulation->GetNode(node_idx);
        c_vector<double, 2> const node_location = node->rGetLocation();
        c_vector<double, 2> const node_to_center = node_location - m_circle_center;
        double const radius = norm(node_to_center);
        if (radius > m_circle_radius)
        {
            c_vector<double, 2> const new_location = m_circle_center + (m_circle_radius - 1e-5) * normalize(node_to_center);
            node->rGetModifiableLocation() = new_location;
        }
    }
}

bool InelasticalCircleBoundaryCondition::VerifyBoundaryCondition()
{
    bool satisfied = true;
    unsigned const num_node = mpCellPopulation->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_node; ++node_idx)
    {
        Node<2>* node = mpCellPopulation->GetNode(node_idx);
        c_vector<double, 2> const node_location = node->rGetLocation();
        double const radius = norm(node_location - m_circle_center);
        if (radius > m_circle_radius)
        {
            satisfied = false;
            break;
        }
    }
    return satisfied;
}

void InelasticalCircleBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
}
