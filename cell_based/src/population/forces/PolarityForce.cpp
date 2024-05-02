#include "PolarityForce.hpp"

PolarityForce::PolarityForce(double self_propelling, std::string polarity_name) : mSelfPropellingParameter(self_propelling),
                                                                                  mPolarityName(std::move(polarity_name)) {}

PolarityForce::~PolarityForce()
{
}

void PolarityForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("MyNagaiHondaForce expects to be used with a VertexBasedCellPopulation");
    }

    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_nodes = p_cell_population->GetNumNodes();
    unsigned const num_elements = p_cell_population->GetNumElements();
    std::vector<c_vector<double, 2> > elem_polarities(num_elements, zero_vector<double>(2));

    for (typename VertexMesh<2, 2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        double theta = 0.0;
        theta = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->GetItem(mPolarityName);
        elem_polarities[elem_idx][0] = cos(theta);
        elem_polarities[elem_idx][1] = sin(theta);
    }

    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        c_vector<double, 2> active_force = zero_vector<double>(2);
        std::set<unsigned> const& relavant_elem_indices = p_cell_population->GetNode(node_idx)->rGetContainingElementIndices();
        for (std::set<unsigned>::const_iterator iter = relavant_elem_indices.begin();
             iter != relavant_elem_indices.end();
             ++iter)
        {
            VertexElement<2, 2>* p_element = p_cell_population->GetElement(*iter);
            double const scalar = mSelfPropellingParameter / static_cast<double>(p_element->GetNumNodes());
            unsigned const elem_idx = p_element->GetIndex();
            active_force += scalar * elem_polarities[elem_idx];
        }
        p_cell_population->GetNode(node_idx)->AddAppliedForceContribution(active_force);
    }
}

void PolarityForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SelfPropellingParameter>" << mSelfPropellingParameter << "</SelfPropellingParameter>\n";
    AbstractForce<2>::OutputForceParameters(rParamsFile);
}