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
#include "ParameterReader.hpp"

void ParamterizedTestPolarityForceWithRandomPolarity(std::string const& path)
{
    try
    {
        FileFinder file_finder{ path, RelativeTo::ChasteSourceRoot };
        ParameterReader reader{};
        reader.parse(file_finder.GetAbsolutePath());

        /*
         *  cell_col: 5
         *  cell_row: 5
         *  sampling_time_multiple: 200
         *  end_time: 20
         *  nagai_honda_deformation_energy_parameter: 100
         *  nagai_honda_membrane_surface_energy_parameter: 10
         *  self_propelling_parameter: 1
         *  target_area_parent_normal_mean: 1
         *  target_area_parent_normal_variance: 1
         *  polarity_delta: 1
         */

        READ_PARAMETER(reader, unsigned, cell_col);
        READ_PARAMETER(reader, unsigned, cell_row);
        READ_PARAMETER(reader, unsigned, sampling_time_multiple);
        READ_PARAMETER(reader, double, end_time);
        READ_PARAMETER(reader, double, nagai_honda_deformation_energy_parameter);
        READ_PARAMETER(reader, double, nagai_honda_membrane_surface_energy_parameter);
        READ_PARAMETER(reader, double, self_propelling_parameter);
        READ_PARAMETER(reader, double, target_area_parent_normal_mean);
        READ_PARAMETER(reader, double, target_area_parent_normal_variance);
        READ_PARAMETER(reader, double, polarity_delta);

        HoneycombVertexMeshGenerator generator{ cell_col, cell_row };
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells{};

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population{ *p_mesh, cells };

        OffLatticeSimulation<2> simulator{ cell_population };
        simulator.SetOutputDirectory(file_finder.GetLeafName());
        simulator.SetSamplingTimestepMultiple(sampling_time_multiple);
        simulator.SetEndTime(end_time);

        MAKE_PTR(MyNagaiHondaForce<2>, p_nagai_honda_force);
        p_nagai_honda_force->SetDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_nagai_honda_force->SetMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        simulator.AddForce(p_nagai_honda_force);

        MAKE_PTR(PolarityForce, p_polarity_force);
        p_polarity_force->SetSelfPropellingParameter(self_propelling_parameter);
        simulator.AddForce(p_polarity_force);

        MAKE_PTR(FixedTargetAreaModifier<2>, p_fixed_target_area_modifier);
        p_fixed_target_area_modifier->SetParentNormalMean(target_area_parent_normal_mean);
        p_fixed_target_area_modifier->SetParentNormalVariance(target_area_parent_normal_variance);
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);

        MAKE_PTR(SimplePolarityModifier<2>, p_simple_polarity_modifier);
        p_simple_polarity_modifier->SetPolarityDelta(polarity_delta);
        simulator.AddSimulationModifier(p_simple_polarity_modifier);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), cell_col * cell_row);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
    catch (Exception& exp)
    {
        TS_FAIL(exp.GetMessage());
    }
}

class TestPolarityForceWithRandomPolarity : public AbstractCellBasedTestSuite
{
public:
    void TestPolarityForceWithRandomPolarity1()
    {
        ParamterizedTestPolarityForceWithRandomPolarity("cell_based/test/pozero/data/PolarityForceWithRandomPolarity1");
    }

    void TestPolarityForceWithRandomPolarity2()
    {
        ParamterizedTestPolarityForceWithRandomPolarity("cell_based/test/pozero/data/PolarityForceWithRandomPolarity2");
    }

    void TestPolarityForceWithRandomPolarity3()
    {
        ParamterizedTestPolarityForceWithRandomPolarity("cell_based/test/pozero/data/PolarityForceWithRandomPolarity3");
    }

    void TestPolarityForceWithRandomPolarity4()
    {
        ParamterizedTestPolarityForceWithRandomPolarity("cell_based/test/pozero/data/PolarityForceWithRandomPolarity4");
    }
};

#endif /*TESTMYNAGAIHONDAFORCEPOLARITYFORCEWITHRANDOMPOLARITY_HPP_*/
