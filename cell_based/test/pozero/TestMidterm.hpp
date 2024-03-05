#ifndef TESTMIDTERM_HPP_
#define TESTMIDTERM_HPP_

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
#include "BoxBoundaryCondition.hpp"
#include "MyNagaiHondaForce.hpp"
#include "FixedTargetAreaModifier.hpp"
#include "PolarityForce.hpp"
#include "PolarityModifier.hpp"
#include "VertexVelocityCorrelationTrackingModifier.hpp"
#include "VertexDisplacementTrackingModifier.hpp"
#include "ParameterReader.hpp"

void ParameterizedTest(std::string const& path)
{
    try
    {
        FileFinder file_finder{ path, RelativeTo::ChasteSourceRoot };
        ParameterReader reader{};
        reader.parse(file_finder.GetAbsolutePath());

        READ_PARAMETER(reader, unsigned, cell_col);
        READ_PARAMETER(reader, unsigned, cell_row);
        READ_PARAMETER(reader, unsigned, sampling_time_multiple);
        READ_PARAMETER(reader, double, end_time);
        READ_PARAMETER(reader, double, nagai_honda_deformation_energy_parameter);
        READ_PARAMETER(reader, double, nagai_honda_membrane_surface_energy_parameter);
        READ_PARAMETER(reader, double, self_propelling_parameter);
        READ_PARAMETER(reader, bool, uniform_target_area);
        READ_PARAMETER(reader, double, initial_area);
        READ_PARAMETER(reader, double, target_area_parent_normal_mean);
        READ_PARAMETER(reader, double, target_area_parent_normal_variance);
        READ_PARAMETER(reader, double, shape_index);
        READ_PARAMETER(reader, double, neightbor_alignment_intensity);
        READ_PARAMETER(reader, double, shape_alignment_intensity);
        READ_PARAMETER(reader, double, protrusion_alignment_intensity);
        READ_PARAMETER(reader, double, correlation_sample_bucket_width);

        HoneycombVertexMeshGenerator generator{ cell_col, cell_row, false, 0.01, 0.001, initial_area };
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

        c_vector<double, 2> aabb_min{};
        c_vector<double, 2> aabb_max{};
        aabb_min[0] = std::numeric_limits<double>::max();
        aabb_min[1] = std::numeric_limits<double>::max();
        aabb_max[0] = std::numeric_limits<double>::min();
        aabb_max[1] = std::numeric_limits<double>::min();
        for (auto iter = p_mesh->GetNodeIteratorBegin(); iter != p_mesh->GetNodeIteratorEnd(); ++iter)
        {
            auto const& position = iter->rGetLocation();
            aabb_min[0] = std::min(aabb_min[0], position[0]);
            aabb_min[1] = std::min(aabb_min[1], position[1]);
            aabb_max[0] = std::max(aabb_max[0], position[0]);
            aabb_max[1] = std::max(aabb_max[1], position[1]);
        }
        auto p_box_boundary_condition = boost::make_shared<BoxBoundaryCondition>(&cell_population, aabb_max[1], aabb_min[1], aabb_min[0], aabb_max[0]);
        simulator.AddCellPopulationBoundaryCondition(p_box_boundary_condition);

        MAKE_PTR(MyNagaiHondaForce<2>, p_nagai_honda_force);
        p_nagai_honda_force->SetDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_nagai_honda_force->SetMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        simulator.AddForce(p_nagai_honda_force);

        MAKE_PTR(PolarityForce, p_polarity_force);
        p_polarity_force->SetSelfPropellingParameter(self_propelling_parameter);
        simulator.AddForce(p_polarity_force);

        boost::shared_ptr<FixedTargetAreaModifier<2> > p_fixed_target_area_modifier{};
        if (uniform_target_area)
        {
            p_fixed_target_area_modifier = boost::make_shared<FixedTargetAreaModifier<2> >(cell_col * cell_row, target_area_parent_normal_mean, shape_index);
        }
        else
        {
            p_fixed_target_area_modifier = boost::make_shared<FixedTargetAreaModifier<2> >(cell_col * cell_row, target_area_parent_normal_mean, target_area_parent_normal_variance, shape_index);
        }
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);

        MAKE_PTR(PolarityModifier, p_polarity_modifier);
        p_polarity_modifier->SetNeighborAlignmentIntensity(neightbor_alignment_intensity);
        p_polarity_modifier->SetShapeAlignmentIntensity(shape_alignment_intensity);
        p_polarity_modifier->SetProtrusionAlignmentIntensity(protrusion_alignment_intensity);
        p_polarity_modifier->SetDt(simulator.GetDt());
        simulator.AddSimulationModifier(p_polarity_modifier);

        auto p_vertex_velocity_correlation_tracker = boost::make_shared<VertexVelocityCorrelationTrackingModifier>(correlation_sample_bucket_width, file_finder.GetLeafName() + "VertexVelocityCorrelation");
        simulator.AddSimulationModifier(p_vertex_velocity_correlation_tracker);

        auto p_vertex_displacement_tracker = boost::make_shared<VertexDisplacementTrackingModifier>(simulator.GetDt(), file_finder.GetLeafName() + "VertexDisplacement");
        simulator.AddSimulationModifier(p_vertex_displacement_tracker);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), cell_col * cell_row);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
    catch (Exception& exp)
    {
        TS_FAIL(exp.GetMessage());
    }
}

class TestMidterm : public AbstractCellBasedTestSuite
{
public:
    void TestMidtermNonUniformTargetArea07() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea07"); }
    // void TestMidtermNonUniformTargetArea08() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea08"); }
    // void TestMidtermNonUniformTargetArea09() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea09"); }
    // void TestMidtermNonUniformTargetArea10() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea10"); }
    // void TestMidtermNonUniformTargetArea11() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea11"); }
    // void TestMidtermNonUniformTargetArea12() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea12"); }
    void TestMidtermNonUniformTargetArea13() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTargetArea13"); }
    // void TestMidtermNonUniformTAVariance01() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTAVariance01"); }
    // void TestMidtermNonUniformTAVariance02() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTAVariance02"); }
    // void TestMidtermNonUniformTAVariance03() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTAVariance03"); }
    // void TestMidtermNonUniformTAVariance04() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTAVariance04"); }
    // void TestMidtermNonUniformTAVariance05() { ParameterizedTest("cell_based/test/pozero/data/MidtermNonUniformTAVariance05"); }
    // void TestMidtermUniformTargetArea07() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea07"); }
    // void TestMidtermUniformTargetArea08() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea08"); }
    // void TestMidtermUniformTargetArea09() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea09"); }
    // void TestMidtermUniformTargetArea10() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea10"); }
    // void TestMidtermUniformTargetArea11() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea11"); }
    // void TestMidtermUniformTargetArea12() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea12"); }
    // void TestMidtermUniformTargetArea13() { ParameterizedTest("cell_based/test/pozero/data/MidtermUniformTargetArea13"); }
};

#endif
