/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTVORONOIVERTEXMESHGENERATOR_HPP_
#define TESTVORONOIVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "VoronoiVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

class TestVoronoiVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestSimpleMesh() throw(Exception)
    {
        // 20 cells wide, 12 high, 4 Lloyd's Relaxation steps, target average element area 1.23
        VoronoiVertexMeshGenerator generator(20, 12, 4, 1.23);
        MutableVertexMesh<2,2>* p_mesh_a = generator.GetMesh();
        MutableVertexMesh<2,2>* p_mesh_b = generator.GetMeshAfterReMesh();

        // Check basic mesh properties are correct
        TS_ASSERT_EQUALS(p_mesh_a->GetNumNodes(), 550u);
        TS_ASSERT_EQUALS(p_mesh_a->GetNumElements(), 240u);

        TS_ASSERT_EQUALS(p_mesh_b->GetNumNodes(), 550u);
        TS_ASSERT_EQUALS(p_mesh_b->GetNumElements(), 240u);

        // Check average cell area is correct
        double average_area = 0.0;
        for (unsigned elem_idx = 0 ; elem_idx < p_mesh_a->GetNumElements() ; elem_idx++)
        {
            average_area += p_mesh_a->GetVolumeOfElement(elem_idx);
        }
        average_area /= double(p_mesh_a->GetNumElements());

        TS_ASSERT_DELTA(average_area, 1.23, 1e-6);
    }

    void TestBoundaryNodes() throw(Exception)
    {
        // 3 cells wide, 2 high, 3 Lloyd's Relaxation steps, target average element area 100.0
        VoronoiVertexMeshGenerator generator(3, 2, 3, 100.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Check basic mesh properties are correct
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 23u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 6u);

        // There should be three non-boundary nodes in this mesh (indices 1, 2, 11)
        TS_ASSERT( !(p_mesh->GetNode(1u)->IsBoundaryNode()) );
        TS_ASSERT( !(p_mesh->GetNode(2u)->IsBoundaryNode()) );
        TS_ASSERT( !(p_mesh->GetNode(11u)->IsBoundaryNode()) );

        // Check that there are exactly 20 boundary nodes
        unsigned num_boundary_nodes = 0;
        for (unsigned node_idx = 0 ; node_idx < 23u ; node_idx++ )
        {
            if (p_mesh->GetNode(node_idx)->IsBoundaryNode())
            {
                num_boundary_nodes++;
            }
        }
        TS_ASSERT_EQUALS(num_boundary_nodes, 20u);
    }

    void TestConstructorExceptions() throw(Exception)
    {
        // Throws because first parameter < 2
        TS_ASSERT_THROWS_THIS(VoronoiVertexMeshGenerator generator(1, 9, 2, 1.23),
                "Need at least 2 by 2 cells");

        // Throws because second parameter < 2
        TS_ASSERT_THROWS_THIS(VoronoiVertexMeshGenerator generator(9, 1, 2, 1.23),
                "Need at least 2 by 2 cells");

        // Throws because fourth parameter <= 0.0
        TS_ASSERT_THROWS_THIS(VoronoiVertexMeshGenerator generator(9, 9, 2, -1.23),
                "Specified target area must be strictly positive");
    }

    void TestValidateSeedLocations() throw(Exception)
    {
        // The instance of the generator class
        VoronoiVertexMeshGenerator generator(3, 3, 0, 1.0);

        // A helper vector
        c_vector<double, 2> one_one;
        one_one[0] = 1.0;
        one_one[1] = 1.0;

        // The sampling multiplier that should be used by the generator
        double sampling_multiplier = 0.5 * double(INT_MAX);

        // Vector of points that will be passed to ValidateSeedLocations()
        std::vector<c_vector<double, 2> > points;

        /*
         * Test two points at 0,0 are moved as expected
         */
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));

        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 0.0, 1e-10);

        /*
         * Test three points at 0,0 are moved as expected
         */
        points.clear();
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));
        points.push_back(zero_vector<double>(2));

        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[2][0], 3.0 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[2][1], 0.0, 1e-10);

        /*
         * Test that periodicity is working correctly
         */
        points.clear();
        points.push_back(one_one);
        points.push_back(one_one);
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / sampling_multiplier, 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.0, 1e-10);

        /*
         * Test that two close non-equal points are moved correctly
         */
        points.clear();
        points.push_back(zero_vector<double>(2));
        points.push_back( (2.0 * DBL_EPSILON) * one_one );
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 0.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);

        /*
         * Test that periodicity for two close non-equal points is working as expected
         */
        points.clear();
        points.push_back(one_one);
        points.push_back((1.0 + 2.0 * DBL_EPSILON) * one_one);
        generator.ValidateSeedLocations(points);

        TS_ASSERT_DELTA(points[0][0], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[0][1], 1.0, 1e-10);
        TS_ASSERT_DELTA(points[1][0], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
        TS_ASSERT_DELTA(points[1][1], 1.5 / (sampling_multiplier * sqrt(2.0)), 1e-10);
    }
};

#endif /*TESTVORONOIVERTEXMESHGENERATOR_HPP_*/
