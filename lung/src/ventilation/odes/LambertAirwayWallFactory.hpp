/*

Copyright (c) 2005-2022, University of Oxford.
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


#ifndef LAMBERTAIRWAYWALLFACTORY_HPP_
#define LAMBERTAIRWAYWALLFACTORY_HPP_

#include "AbstractAirwayWallFactory.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "AirwayTreeWalker.hpp"
#include "LambertAirwayWall.hpp"

/**
 * A factory to ease creating Lambert based airway wall models.
 *
 * This creates Lambert models based on Horsfield or Strahler order.
 * Generation data from Lambert's 1982 paper is interpolated onto the
 * corresponding order.
 */
class LambertAirwayWallFactory : public AbstractAirwayWallFactory
{
public:

    /**
     * Default constructor.
     *
     * By default airway walls models will be defined by interpolating Lambert's parameters
     * onto element Horsfield order.
     *
     * @param useStrahlerOrder Choose to use Strahler order instead of Horsfield order
     */
    LambertAirwayWallFactory(bool useStrahlerOrder = false);

    /**
     * Destructor
     */
    virtual ~LambertAirwayWallFactory();

    /**
     * @return a newly airway wall object for the given node.
     *
     * @param pNode  Pointer to node object.
     * @return An airway wall object parameterised for this airway
     */
    virtual LambertAirwayWall* CreateAirwayWallForElement(Element<1,3>* pElement);

    double GetAlpha0ForGeneration(unsigned generation);
    double GetAlpha0PrimeForGeneration(unsigned generation);
    double GetN1ForGeneration(unsigned generation);
    double GetN2ForGeneration(unsigned generation);

    /**
     * @param time The current time in seconds
     * @param pNode Pointer to node object.
     * @return The pleural pressure at the given node at the given time
     */
    virtual double GetPleuralPressureForAirway(double time, Element<1,3>* pElement);

    /**
     * @param pMesh  the mesh for which to create acinar units.
     */
    virtual void SetMesh(AbstractTetrahedralMesh<1,3>* pMesh);

    /** alpha_0 values by generation dimensionless */
    static const double mAlpha0[17];

    /** alpha_0_prime values by generation in Pascals^-1 */
    static const double mAlpha0Prime[17];

    /** n1 values by generation dimensionless */
    static const double mN1[17];

    /** n2 values by generation dimensionless  */
    static const double mN2[17];

    /** The number of generations in the Lambert data. Nb generation number starts at 0! */
    static const double mMaxGeneration;

private:
    /** Walker to determine order of airway elements */
    AirwayTreeWalker* mpWalker;

    /** By default use Horsfield order, this flag allows use of Strahler order */
    bool mUseStrahlerOrder;

    /** The largest order in the airway tree*/
    double mMaxOrder;
};

#endif /*LAMBERTAIRWAYWALLFACTORY_HPP_*/
