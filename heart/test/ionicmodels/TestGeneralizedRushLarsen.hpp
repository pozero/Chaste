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

#ifndef TESTGENERALIZEDRUSHLARSEN_HPP_
#define TESTGENERALIZEDRUSHLARSEN_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp" // Needed to avoid memory error on Boost 1.34

#include "LuoRudy1991.hpp"
#include "LuoRudy1991Opt.hpp"
#include "AbstractGeneralizedRushLarsenCardiacCell.hpp" // Needed for chaste_libs=0 build
#include "ZeroStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "GRL1IvpOdeSolver.hpp"
#include "GRL2IvpOdeSolver.hpp"

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"
#include "CellMLToSharedLibraryConverter.hpp"
#include "DynamicCellModelLoader.hpp"
#include "Warnings.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
Megan E. Marsh, Raymond J. Spiteri
Numerical Simulation Laboratory
University of Saskatchewan
December 2011
Partial support provided by research grants from the National
Science and Engineering Research Council (NSERC) of Canada
and the MITACS/Mprime Canadian Network of Centres of Excellence.
*/


/**
 * This test is based on code TestRushLarsen and tests the GRL1 and GRL2
 * autogenerated Luo--Rudy code
 */
class TestGeneralizedRushLarsen : public CxxTest::TestSuite
{
    AbstractCardiacCell* mpGeneralizedRushLarsenCell;
//    AbstractCardiacCell* mpGeneralizedRushLarsenCellOpt;

    void GenerateCells()
    {
        // Do the conversions preserving generated sources
        CellMLToSharedLibraryConverter converter(true);
        std::string dirname = "TestGeneralizedRushLarsen";
        std::string model = "LuoRudy1991";
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());

        std::vector<std::string> args;
        args.push_back("--grl1");

        { // No opt
            // Copy CellML file into output dir
            OutputFileHandler handler(dirname + "/normal");
            FileFinder cellml_file("heart/src/odes/cellml/" + model + ".cellml", RelativeTo::ChasteSourceRoot);
            FileFinder copied_file = handler.CopyFileTo(cellml_file);

            // Create options file & convert
            converter.SetOptions(args);
            DynamicCellModelLoaderPtr p_loader = converter.Convert(copied_file);
            mpGeneralizedRushLarsenCell = dynamic_cast<AbstractCardiacCell*>(p_loader->CreateCell(p_solver, p_stimulus));
        }
    }

    void GenerateCells2()
    {
        // Do the conversions preserving generated sources
        CellMLToSharedLibraryConverter converter(true);
        std::string dirname = "TestGeneralizedRushLarsen2";
        std::string model = "LuoRudy1991";
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());

        std::vector<std::string> args;
        args.push_back("--grl2");

        { // No opt
            // Copy CellML file into output dir
            OutputFileHandler handler(dirname + "/normal");
            FileFinder cellml_file("heart/src/odes/cellml/" + model + ".cellml", RelativeTo::ChasteSourceRoot);
            FileFinder copied_file = handler.CopyFileTo(cellml_file);

            // Create options file & convert
            converter.SetOptions(args);
            DynamicCellModelLoaderPtr p_loader = converter.Convert(copied_file);
            mpGeneralizedRushLarsenCell = dynamic_cast<AbstractCardiacCell*>(p_loader->CreateCell(p_solver, p_stimulus));
        }
    }

public:
    void TestLuoRudyGeneralizedRushLarsenMethod()
    {
        EXIT_IF_PARALLEL;
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        GenerateCells();

        // Check the models really use Rush-Larsen
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        // Some coverage
        AbstractGeneralizedRushLarsenCardiacCell* p_grl_cell = dynamic_cast<AbstractGeneralizedRushLarsenCardiacCell*>(mpGeneralizedRushLarsenCell);
        TS_ASSERT(p_grl_cell);
        TS_ASSERT(!p_grl_cell->HasAnalyticJacobian());
        TS_ASSERT(!p_grl_cell->GetUseAnalyticJacobian());
        TS_ASSERT_THROWS_THIS(p_grl_cell->ForceUseOfNumericalJacobian(false),
                              "Using analytic Jacobian terms for generalised Rush-Larsen is not yet supported.");
        TS_ASSERT_THROWS_NOTHING(p_grl_cell->ForceUseOfNumericalJacobian(true));

        // Normal Luo-Rudy for comparison
        boost::shared_ptr<HeunIvpOdeSolver> p_heun_solver(new HeunIvpOdeSolver());
        CellLuoRudy1991FromCellML reference_model(p_heun_solver, mpGeneralizedRushLarsenCell->GetStimulusFunction());

        // Check GetIIonic is identical
        TS_ASSERT_DELTA(mpGeneralizedRushLarsenCell->GetIIonic(), reference_model.GetIIonic(), 1e-12);

        // Test non-stimulated cell (using ComputeExceptVoltage)
        mpGeneralizedRushLarsenCell->ComputeExceptVoltage(0.0, 1.0);
        reference_model.ComputeExceptVoltage(0.0, 1.0);

        TS_ASSERT_EQUALS(mpGeneralizedRushLarsenCell->GetNumberOfStateVariables(),
                         reference_model.GetNumberOfStateVariables());
        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(mpGeneralizedRushLarsenCell->rGetStateVariables()[i],
                            reference_model.rGetStateVariables()[i], 1e-6);
        }

        // Test stimulated cell (using Compute)
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-25.5, 1.99, 0.0));
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        mpGeneralizedRushLarsenCell->SetStimulusFunction(p_stimulus);
        OdeSolution solutions_GRL1 = mpGeneralizedRushLarsenCell->Compute(0.0, 1.0, 0.01);
        TS_ASSERT_EQUALS(solutions_GRL1.GetNumberOfTimeSteps(), 100u);

        reference_model.ResetToInitialConditions();
        reference_model.SetStimulusFunction(p_stimulus);
        OdeSolution solutions_ref = reference_model.Compute(0.0, 1.0, 0.01);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(solutions_GRL1.rGetSolutions().back()[i],
                            solutions_ref.rGetSolutions().back()[i], 1e-2);
        }

        //Compare LuoRudy solution with general GRL1 solver solution (should match)
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        boost::shared_ptr<ZeroStimulus> p_stimulus_zero(new ZeroStimulus());
        mpGeneralizedRushLarsenCell->SetStimulusFunction(p_stimulus_zero);
        mpGeneralizedRushLarsenCell->SetTimestep(1e-3);
        mpGeneralizedRushLarsenCell->SetVoltage(-30);
        OdeSolution solutions_GRL1_stimulated_cell_order1 = mpGeneralizedRushLarsenCell->Compute(0.0, 1e-3);

        // Compare with SolveAndUpdateState for coverage
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        mpGeneralizedRushLarsenCell->SetVoltage(-30);
        mpGeneralizedRushLarsenCell->SolveAndUpdateState(0.0, 1e-3);
        for (unsigned i=0; i<mpGeneralizedRushLarsenCell->GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(mpGeneralizedRushLarsenCell->rGetStateVariables()[i],
                            solutions_GRL1_stimulated_cell_order1.rGetSolutions().back()[i], 1e-12);
        }

        // Compare with general GRL1 method
        boost::shared_ptr<GRL1IvpOdeSolver> p_grl1_solver(new GRL1IvpOdeSolver());
        CellLuoRudy1991FromCellML reference_model_grl1(p_grl1_solver, mpGeneralizedRushLarsenCell->GetStimulusFunction());
        reference_model_grl1.ResetToInitialConditions();
        reference_model_grl1.SetStimulusFunction(p_stimulus_zero);
        reference_model_grl1.SetTimestep(1e-3);
        reference_model_grl1.SetVoltage(-30);
        OdeSolution ref_solution_grl1 = reference_model_grl1.Compute(0.0, 1e-3);

        for (unsigned i=0; i<reference_model_grl1.GetNumberOfStateVariables(); i++)
        {
            double error1 = solutions_GRL1_stimulated_cell_order1.rGetSolutions().back()[i] - ref_solution_grl1.rGetSolutions().back()[i];
            TS_ASSERT_DELTA(fabs(error1),0,5e-7);
        }

        // Test order of convergence (first-order method)

        // Get two solutions with halved stepsize
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        reference_model_grl1.SetStimulusFunction(p_stimulus_zero);
        reference_model_grl1.SetVoltage(-30);
        mpGeneralizedRushLarsenCell->SetTimestep(0.002);
        solutions_GRL1_stimulated_cell_order1 = mpGeneralizedRushLarsenCell->Compute(0.0, 1.0);

        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        reference_model_grl1.SetStimulusFunction(p_stimulus_zero);
        reference_model_grl1.SetVoltage(-30);
        mpGeneralizedRushLarsenCell->SetTimestep(0.001);
        OdeSolution solutions_GRL1_stimulated_cell_order2 = mpGeneralizedRushLarsenCell->Compute(0.0, 1.0);

        // Get accurate solution to be used as reference solution
        //boost::shared_ptr<HeunIvpOdeSolver> p_heun_solver(new HeunIvpOdeSolver());
        reference_model.SetStimulusFunction(p_stimulus_zero);
        reference_model.SetVoltage(-30);
        reference_model.ResetToInitialConditions();
        reference_model.SetTimestep(1e-6);
        OdeSolution ref_solution = reference_model.Compute(0.0, 1.0);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            double error1 = solutions_GRL1_stimulated_cell_order1.rGetSolutions().back()[i] - ref_solution.rGetSolutions().back()[i];
            double error2 = solutions_GRL1_stimulated_cell_order2.rGetSolutions().back()[i] - ref_solution.rGetSolutions().back()[i];
            TS_ASSERT_DELTA(error1/error2, 2, 0.1);
        }

        // Free memory
        delete mpGeneralizedRushLarsenCell;
//        delete mpGeneralizedRushLarsenCellOpt;
    }

    void TestLuoRudyGeneralizedRushLarsenMethod2()
    {
        EXIT_IF_PARALLEL;
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        GenerateCells2();

        // Check the models really use Rush-Larsen
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        // Normal Luo-Rudy for comparison
        boost::shared_ptr<HeunIvpOdeSolver> p_heun_solver(new HeunIvpOdeSolver());
        CellLuoRudy1991FromCellML reference_model(p_heun_solver, mpGeneralizedRushLarsenCell->GetStimulusFunction());

        // Check GetIIonic is identical
        TS_ASSERT_DELTA(mpGeneralizedRushLarsenCell->GetIIonic(), reference_model.GetIIonic(), 1e-12);

        // Test non-stimulated cell (using ComputeExceptVoltage)
        mpGeneralizedRushLarsenCell->ComputeExceptVoltage(0.0, 1.0);
        reference_model.ComputeExceptVoltage(0.0, 1.0);

        TS_ASSERT_EQUALS(mpGeneralizedRushLarsenCell->GetNumberOfStateVariables(),
                         reference_model.GetNumberOfStateVariables());
        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(mpGeneralizedRushLarsenCell->rGetStateVariables()[i],
                            reference_model.rGetStateVariables()[i], 1e-6);
        }
        // Test stimulated cell (using Compute)
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(-25.5, 1.99, 0.0));
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        mpGeneralizedRushLarsenCell->SetStimulusFunction(p_stimulus);
        OdeSolution solutions_GRL2 = mpGeneralizedRushLarsenCell->Compute(0.0, 1.0, 0.01);
        TS_ASSERT_EQUALS(solutions_GRL2.GetNumberOfTimeSteps(), 100u);

        reference_model.ResetToInitialConditions();
        reference_model.SetStimulusFunction(p_stimulus);
        OdeSolution solutions_ref = reference_model.Compute(0.0, 1.0, 0.01);

        for (unsigned i=0; i<reference_model.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(solutions_GRL2.rGetSolutions().back()[i],
                            solutions_ref.rGetSolutions().back()[i], 1e-2);
        }

        //Compare LuoRudy solution with general GRL2 solver solution (should match)
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        boost::shared_ptr<ZeroStimulus> p_stimulus_zero(new ZeroStimulus());
        mpGeneralizedRushLarsenCell->ResetToInitialConditions();
        mpGeneralizedRushLarsenCell->SetStimulusFunction(p_stimulus_zero);
        mpGeneralizedRushLarsenCell->SetTimestep(1e-3);
        mpGeneralizedRushLarsenCell->SetVoltage(-30);
        OdeSolution solutions_GRL2_stimulated_cell_order2 = mpGeneralizedRushLarsenCell->Compute(0.0, 1e-3);

        boost::shared_ptr<GRL2IvpOdeSolver> p_grl2_solver(new GRL2IvpOdeSolver());
        CellLuoRudy1991FromCellML reference_model_grl2(p_grl2_solver, mpGeneralizedRushLarsenCell->GetStimulusFunction());

        // Compare with general GRL2 method
        reference_model_grl2.ResetToInitialConditions();
        reference_model_grl2.SetStimulusFunction(p_stimulus_zero);
        reference_model_grl2.SetTimestep(1e-3);
        reference_model_grl2.SetVoltage(-30);
        OdeSolution ref_solution_grl2 = reference_model_grl2.Compute(0.0, 1e-3);

        for (unsigned i=0; i<reference_model_grl2.GetNumberOfStateVariables(); i++)
        {
            double error1 = solutions_GRL2_stimulated_cell_order2.rGetSolutions().back()[i] - ref_solution_grl2.rGetSolutions().back()[i];
            TS_ASSERT_DELTA(fabs(error1),0,5e-7);
        }

        // Free memory
        delete mpGeneralizedRushLarsenCell;
//        delete mpGeneralizedRushLarsenCellOpt;
    }
};

#endif // TESTGENERALIZEDRUSHLARSEN_HPP_
