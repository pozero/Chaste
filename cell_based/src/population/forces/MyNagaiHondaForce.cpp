#include "MyNagaiHondaForce.hpp"

template <unsigned DIM>
MyNagaiHondaForce<DIM>::MyNagaiHondaForce() : AbstractForce<DIM>(),
                                              mNagaiHondaDeformationEnergyParameter(100.0),
                                              mNagaiHondaMembraneSurfaceEnergyParameter(10.0)
{
}

template <unsigned DIM>
MyNagaiHondaForce<DIM>::~MyNagaiHondaForce()
{
}

template <unsigned DIM>
void MyNagaiHondaForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("MyNagaiHondaForce expects to be used with a VertexBasedCellPopulation");
    }

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned const num_nodes = p_cell_population->GetNumNodes();
    unsigned const num_elements = p_cell_population->GetNumElements();
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    std::vector<double> target_perimeters(num_elements);

    for (typename VertexMesh<DIM, DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned const elem_idx = elem_iter->GetIndex();
        element_areas[elem_idx] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_idx);
        element_perimeters[elem_idx] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_idx);
        try
        {
            target_areas[elem_idx] = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("MyNagaiHondaForce need target area to proceed");
        }
        try
        {
            target_perimeters[elem_idx] = p_cell_population->GetCellUsingLocationIndex(elem_idx)->GetCellData()->GetItem("target perimeter");
        }
        catch (Exception&)
        {
            EXCEPTION("MyNagaiHondaForce need target perimeter to proceed");
        }
    }

    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> surface_expansion_contribution = zero_vector<double>(DIM);
        std::set<unsigned> const& relavant_elem_indices = p_cell_population->GetNode(node_idx)->rGetContainingElementIndices();
        for (std::set<unsigned>::const_iterator iter = relavant_elem_indices.begin();
             iter != relavant_elem_indices.end();
             ++iter)
        {
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned const elem_idx = p_element->GetIndex();
            unsigned const local_idx = p_element->GetNodeLocalIndex(node_idx);

            c_vector<double, DIM> elem_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_idx);
            deformation_contribution -= mNagaiHondaDeformationEnergyParameter * (element_areas[elem_idx] - target_areas[elem_idx]) * elem_area_gradient;

            c_vector<double, DIM> elem_perimeter_gradient = p_cell_population->rGetMesh().GetPerimeterGradientOfElementAtNode(p_element, local_idx);
            membrane_surface_tension_contribution -= mNagaiHondaMembraneSurfaceEnergyParameter * (element_perimeters[elem_idx] - target_perimeters[elem_idx]) * elem_perimeter_gradient;

            surface_expansion_contribution += mCellExpansionTerm * elem_area_gradient;
        }
        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + surface_expansion_contribution;
        p_cell_population->GetNode(node_idx)->AddAppliedForceContribution(force_on_node);
    }
}

template <unsigned DIM>
double MyNagaiHondaForce<DIM>::GetDeformationEnergyParameter() const
{
    return mNagaiHondaDeformationEnergyParameter;
}

template <unsigned DIM>
double MyNagaiHondaForce<DIM>::GetMembraneSurfaceEnergyParameter() const
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template <unsigned DIM>
void MyNagaiHondaForce<DIM>::SetDeformationEnergyParameter(double newParam)
{
    mNagaiHondaDeformationEnergyParameter = newParam;
}

template <unsigned DIM>
void MyNagaiHondaForce<DIM>::SetMembraneSurfaceEnergyParameter(double newParam)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = newParam;
}

template <unsigned DIM>
void MyNagaiHondaForce<DIM>::SetCellExpansionTerm(double newParam)
{
    mCellExpansionTerm = newParam;
}

template <unsigned DIM>
void MyNagaiHondaForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template class MyNagaiHondaForce<1>;
template class MyNagaiHondaForce<2>;
template class MyNagaiHondaForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNagaiHondaForce);
