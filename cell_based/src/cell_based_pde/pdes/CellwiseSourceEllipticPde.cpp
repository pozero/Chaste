/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellwiseSourceEllipticPde.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CellwiseSourceEllipticPde<DIM>::CellwiseSourceEllipticPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double coefficient)
    : mrCellPopulation(rCellPopulation),
      mCoefficient(coefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& CellwiseSourceEllipticPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::GetCoefficient() const
{
    return mCoefficient;
}

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return 0.0;
}

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}

///\todo Consider how to reduce the number of dynamic_casts, here and elsewhere (#2687)
template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    double coefficient = 0.0;
    unsigned tet_node_index = rNode.GetIndex();
    bool is_cell_apoptotic = false;

    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)) != NULL ||
        dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)) != NULL)  // Intel compiler wants the "!= NULL"
    {
        if (this->mrCellPopulation.IsCellAttachedToLocationIndex(tet_node_index))
        {
            /*
             * For a centre-based cell population, this node of the tetrahedral finite element mesh
             * is the same as the node in the population and associated with a cell; for a cellular
             * Potts cell population, this node corresponds to the Potts element associated with
             * a cell.
             */
            is_cell_apoptotic = this->mrCellPopulation.GetCellUsingLocationIndex(tet_node_index)->template HasCellProperty<ApoptoticCellProperty>();
        }
        else
        {
            // There is no cell associated with this node
            return 0.0;
        }
    }
    else if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)) != NULL)  // Intel compiler wants the "!= NULL"
    {
        VertexBasedCellPopulation<DIM>* static_cast_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));

        if (rNode.GetIndex() < static_cast_cell_population->GetNumNodes())
        {
            std::set<unsigned> containing_element_indices = static_cast_cell_population->GetNode(tet_node_index)->rGetContainingElementIndices();

            for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
                 iter != containing_element_indices.end();
                 iter++)
            {
                if (static_cast_cell_population->GetCellUsingLocationIndex(*iter)->template HasCellProperty<ApoptoticCellProperty>() )
                {
                    is_cell_apoptotic = true;
                    break;
                }
            }
        }
        else
        {
            /*
             * This node of the tetrahedral finite element mesh is in the centre of the element of the
             * vertex based cell population, so we can use an offset to compute which cell to interrogate.
             */
            is_cell_apoptotic = this->mrCellPopulation.GetCellUsingLocationIndex(rNode.GetIndex()-static_cast_cell_population->GetNumNodes())->template HasCellProperty<ApoptoticCellProperty>();
        }
    }
    else if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)) != NULL)  // Intel compiler wants the "!= NULL"
    {
        /*
         * For a CA-based cell population, the index of this node in the tetrahedral finite element mesh
         * corresponds to the 'position' of the cell to interrogate in the vector of cells.
         */
        typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();

        assert(tet_node_index < this->mrCellPopulation.GetNumRealCells());
        for (unsigned i=0; i<tet_node_index; i++)
        {
            ++cell_iter;
        }
        is_cell_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();
    }
    else
    {
        NEVER_REACHED;
    }

    if (!is_cell_apoptotic)
    {
        coefficient = mCoefficient;
    }

    return coefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> CellwiseSourceEllipticPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}

// Explicit instantiation
template class CellwiseSourceEllipticPde<1>;
template class CellwiseSourceEllipticPde<2>;
template class CellwiseSourceEllipticPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellwiseSourceEllipticPde)
