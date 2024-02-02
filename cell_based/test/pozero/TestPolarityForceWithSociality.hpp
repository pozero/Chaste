#ifndef TESTPOLARITYFORCEWITHSOCIALITY_HPP_
#define TESTPOLARITYFORCEWITHSOCIALITY_HPP_

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
#include "InelasticalCircleBoundaryCondition.hpp"
#include "MyNagaiHondaForce.hpp"
#include "FixedTargetAreaModifier.hpp"
#include "PolarityForce.hpp"
#include "PolarityModifier.hpp"
#include "ParameterReader.hpp"

void ParameterizedTestPolarityForceWithSociality(std::string const& path)
{
    try
    {
        FileFinder file_finder{ path, RelativeTo::ChasteSourceRoot };
        ParameterReader reader{};
        reader.parse(file_finder.GetAbsolutePath());

        /*
            cell_col: 6
            cell_row: 6
            sampling_time_multiple: 50
            end_time: 20
            nagai_honda_deformation_energy_parameter: 10
            nagai_honda_membrane_surface_energy_parameter: 2
            self_propelling_parameter: 1
            target_area_parent_normal_mean: 1
            target_area_parent_normal_variance: 0.2
            shape_index: 1.0
            neightbor_alignment_intensity: 0.5
            shape_alignment_intensity: 1
            protrusion_alignment_intensity: 35
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
        READ_PARAMETER(reader, double, shape_index);
        READ_PARAMETER(reader, double, neightbor_alignment_intensity);
        READ_PARAMETER(reader, double, shape_alignment_intensity);
        READ_PARAMETER(reader, double, protrusion_alignment_intensity);

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

        c_vector<double, 2> boundary_center = zero_vector<double>(2);
        double biggest_width = std::numeric_limits<double>::min();
        biggest_width = std::max(cell_population.GetWidth(0), biggest_width);
        biggest_width = std::max(cell_population.GetWidth(1), biggest_width);
        double const boundary_radius = (sqrt(2.0) / 2.0) * biggest_width;
        boundary_center[0] = boundary_radius * 0.5 + 1.0;
        boundary_center[1] = boundary_radius * 0.5 + 1.0;
        boost::shared_ptr<InelasticalCircleBoundaryCondition> p_inelastical_circle_boundary_condition{
            new InelasticalCircleBoundaryCondition(&cell_population, boundary_center, boundary_radius)
        };
        simulator.AddCellPopulationBoundaryCondition(p_inelastical_circle_boundary_condition);

        MAKE_PTR(MyNagaiHondaForce<2>, p_nagai_honda_force);
        p_nagai_honda_force->SetDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_nagai_honda_force->SetMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        simulator.AddForce(p_nagai_honda_force);

        MAKE_PTR(PolarityForce, p_polarity_force);
        p_polarity_force->SetSelfPropellingParameter(self_propelling_parameter);
        simulator.AddForce(p_polarity_force);

        boost::shared_ptr<FixedTargetAreaModifier<2> > p_fixed_target_area_modifier{ new FixedTargetAreaModifier<2>{ cell_col * cell_row, target_area_parent_normal_mean, target_area_parent_normal_variance, shape_index } };
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);

        MAKE_PTR(PolarityModifier, p_polarity_modifier);
        p_polarity_modifier->SetNeighborAlignmentIntensity(neightbor_alignment_intensity);
        p_polarity_modifier->SetShapeAlignmentIntensity(shape_alignment_intensity);
        p_polarity_modifier->SetProtrusionAlignmentIntensity(protrusion_alignment_intensity);
        p_polarity_modifier->SetDt(simulator.GetDt());
        simulator.AddSimulationModifier(p_polarity_modifier);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), cell_col * cell_row);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
    catch (Exception& exp)
    {
        TS_FAIL(exp.GetMessage());
    }
}

class TestPolarityForceWithSociality : public AbstractCellBasedTestSuite
{
public:
    void TestPolarityForceWithSociality1() { ParameterizedTestPolarityForceWithSociality("cell_based/test/pozero/data/PolarityForceWithSociality1"); }
    void TestPolarityForceWithSociality2() { ParameterizedTestPolarityForceWithSociality("cell_based/test/pozero/data/PolarityForceWithSociality2"); }
    void TestPolarityForceWithSociality3() { ParameterizedTestPolarityForceWithSociality("cell_based/test/pozero/data/PolarityForceWithSociality3"); }
};

#endif
