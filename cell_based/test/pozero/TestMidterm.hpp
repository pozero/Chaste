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

void ParameterizedTest(std::string const& base, std::string const& config)
{
    try
    {
        ParameterReader reader{};
        {
            FileFinder file_finder{ base, RelativeTo::ChasteSourceRoot };
            reader.parse(file_finder.GetAbsolutePath());
        }
        FileFinder file_finder{ config, RelativeTo::ChasteSourceRoot };
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
        READ_PARAMETER(reader, double, target_area_min);
        READ_PARAMETER(reader, double, target_area_max);
        READ_PARAMETER(reader, double, shape_index);
        READ_PARAMETER(reader, double, neightbor_alignment_intensity);
        READ_PARAMETER(reader, double, shape_alignment_intensity);
        READ_PARAMETER(reader, double, protrusion_alignment_intensity);
        READ_PARAMETER(reader, unsigned, correlation_sample_bucket_count);

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
        double const box_diameter = norm_2(aabb_max - aabb_min);
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
            p_fixed_target_area_modifier = boost::make_shared<FixedTargetAreaModifier<2> >(cell_col * cell_row, target_area_parent_normal_mean, target_area_parent_normal_variance, target_area_min, target_area_max, shape_index);
        }
        simulator.AddSimulationModifier(p_fixed_target_area_modifier);

        MAKE_PTR(PolarityModifier, p_polarity_modifier);
        p_polarity_modifier->SetNeighborAlignmentIntensity(neightbor_alignment_intensity);
        p_polarity_modifier->SetShapeAlignmentIntensity(shape_alignment_intensity);
        p_polarity_modifier->SetProtrusionAlignmentIntensity(protrusion_alignment_intensity);
        p_polarity_modifier->SetDt(simulator.GetDt());
        simulator.AddSimulationModifier(p_polarity_modifier);

        auto p_vertex_velocity_correlation_tracker = boost::make_shared<VertexVelocityCorrelationTrackingModifier>(box_diameter, correlation_sample_bucket_count, file_finder.GetLeafName() + "VertexVelocityCorrelation");
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
    void TestBase() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermBase"); }

    // void TestMidtermNonUniformTargetArea17() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea17"); }
    // void TestMidtermNonUniformTargetArea18() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea18"); }
    // void TestMidtermNonUniformTargetArea19() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea19"); }
    // void TestMidtermNonUniformTargetArea20() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea20"); }
    // void TestMidtermNonUniformTargetArea21() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea21"); }
    // void TestMidtermNonUniformTargetArea22() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea22"); }
    // void TestMidtermNonUniformTargetArea23() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTargetArea23"); }

    // void TestMidtermNonUniformTAVariance08() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTAVariance08"); }
    // void TestMidtermNonUniformTAVariance09() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTAVariance09"); }
    // void TestMidtermNonUniformTAVariance10() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTAVariance10"); }
    // void TestMidtermNonUniformTAVariance11() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTAVariance11"); }
    // void TestMidtermNonUniformTAVariance12() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermNonUniformTAVariance12"); }

    // void TestMidtermUniformTargetArea17() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea17"); }
    // void TestMidtermUniformTargetArea18() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea18"); }
    // void TestMidtermUniformTargetArea19() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea19"); }
    // void TestMidtermUniformTargetArea20() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea20"); }
    // void TestMidtermUniformTargetArea21() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea21"); }
    // void TestMidtermUniformTargetArea22() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea22"); }
    // void TestMidtermUniformTargetArea23() { ParameterizedTest("cell_based/test/pozero/data/MidtermBase", "cell_based/test/pozero/data/MidtermUniformTargetArea23"); }
};

#endif
