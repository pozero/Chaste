/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifndef TESTCODEGENLONG_HPP_
#define TESTCODEGENLONG_HPP_

#include "CodegenLongHelperTestSuite.hpp"

#include <boost/foreach.hpp>
#include <vector>

#include "HeartConfig.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test chaste_codegen functionality to generate Normal cells,
 * by dynamically loading (and hence converting) a wide range of cell models.
 */
class TestCodegenLongNormal : public CodegenLongHelperTestSuite
{
public:
    void TestNormalCells()
    {
        std::cout << "Search for 'Failure', ': ***', 'Error', or 'Failed' to find problems." << std::endl;

        std::string dirname("TestCodegenLongNormal");
        std::vector<std::string> args; 

        // Models that need a very small dt
        std::vector<std::string> small_dt_models = spectail_streatment_models(models, {"li_mouse_2010"});

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname, models, args);

        // See Cooper Spiteri Mirams paper table 2
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0001953125, 0.1, 1.0);
        RunTests(dirname + "-small-dt", small_dt_models, args);
    }
};

#endif // TESTCODEGENLONGNORMAL_HPP_
