#ifndef TESTVORONOIGENERATOR_HPP_
#define TESTVORONOIGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "MyVoronoiVertexMeshGenerator.hpp"
#include "Cell.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "RandomNumberGenerator.hpp"

class TestVoronoiGenerator : public AbstractCellBasedTestSuite
{
public:
    void Test1()
    {
        try
        {
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            std::vector<Vector2> sites(5 * 5, Vector2{});
            for (Vector2& site : sites)
            {
                site.x = p_gen->ranf();
                site.y = p_gen->ranf();
            }
            Box bounding_box{ -0.05, -0.05, 1.05, 1.05 };
            MyVoronoiVertexMeshGenerator mesh_generator{ sites, bounding_box };
            MutableVertexMesh<2, 2>* p_mesh = mesh_generator.GetMesh();
            std::vector<CellPtr> cells{};

            MAKE_PTR(TransitCellProliferativeType, p_transit_type);
            CellsGenerator<NoCellCycleModel, 2> cells_generator;

            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

            VertexBasedCellPopulation<2> cell_population{ *p_mesh, cells };
            OffLatticeSimulation<2> simulator{ cell_population };
            simulator.SetOutputDirectory("MyVoronoiGenerator");
            simulator.SetSamplingTimestepMultiple(200);
            simulator.SetEndTime(20.0);
            simulator.Solve();

            TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25);
        }
        catch (Exception& exp)
        {
            TS_FAIL(exp.GetMessage());
        }
    }
};

#endif
