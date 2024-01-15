#ifndef TESTMYNAGAIHONDAFORCEPOLARITYFORCEWITHRANDOMPOLARITY_HPP_
#define TESTMYNAGAIHONDAFORCEPOLARITYFORCEWITHRANDOMPOLARITY_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "Cell.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "MyNagaiHondaForce.hpp"
#include "FixedTargetAreaModifier.hpp"
#include "PolarityForce.hpp"
#include "SimplePolarityModifier.hpp"

class TestMyNagaiHondaForcePolarityForceWithRandomPolarity : public AbstractCellBasedTestSuite
{
public:
    void TestMyNagaiHondaForce()
    {
        HoneycombVertexMeshGenerator generator{ 5, 5 };
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells{};

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population{ *p_mesh, cells };

        OffLatticeSimulation<2> simulator{ cell_population };
        simulator.SetOutputDirectory("MyNagiHondaForceOnly");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);
        MAKE_PTR(MyNagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);
        MAKE_PTR(FixedTargetAreaModifier<2>, p_fixed_target_area_modifier);
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    void TestMyNagaiHondaForcePolarityForce()
    {
        HoneycombVertexMeshGenerator generator{ 5, 5 };
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells{};

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population{ *p_mesh, cells };

        OffLatticeSimulation<2> simulator{ cell_population };
        simulator.SetOutputDirectory("MyNagiHondaForcePolarityForce");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        MAKE_PTR(MyNagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);
        MAKE_PTR(PolarityForce, p_polarity_force);
        simulator.AddForce(p_polarity_force);

        MAKE_PTR(FixedTargetAreaModifier<2>, p_fixed_target_area_modifier);
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);
        MAKE_PTR(SimplePolarityModifier<2>, p_simple_polarity_modifier);
        simulator.AddSimulationModifier(p_simple_polarity_modifier);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }
};

#endif /*TESTMYNAGAIHONDAFORCEPOLARITYFORCEWITHRANDOMPOLARITY_HPP_*/
