#ifndef TESTFINAL_HPP_
#define TESTFINAL_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "BoxBoundaryCondition.hpp"
#include "NagaiHondaForce.hpp"
#include "PolarityForce.hpp"
#include "PolarityModifier.hpp"
#include "SimplePolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "MyVoronoiVertexMeshGenerator.hpp"
#include "Cell.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "MyNagaiHondaForce.hpp"
#include "FixedTargetAreaModifier.hpp"
#include "ParameterReader.hpp"
#include "VertexDisplacementTrackingModifier.hpp"
#include "VertexVelocityCorrelationTrackingModifier.hpp"

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
        READ_PARAMETER(reader, unsigned, random_seed);
        READ_PARAMETER(reader, unsigned, sampling_time_multiple);
        READ_PARAMETER(reader, double, end_time);
        READ_PARAMETER(reader, double, nagai_honda_deformation_energy_parameter);
        READ_PARAMETER(reader, double, nagai_honda_membrane_surface_energy_parameter);
        READ_PARAMETER(reader, double, nagai_honda_cell_expansion_parameter);
        READ_PARAMETER(reader, double, align_propelling_parameter);
        READ_PARAMETER(reader, double, random_propelling_parameter);
        READ_PARAMETER(reader, unsigned, cell_count);
        READ_PARAMETER(reader, unsigned, uniform_taret_area);
        READ_PARAMETER(reader, double, target_area_parent_normal_mean);
        READ_PARAMETER(reader, double, target_area_parent_normal_variance);
        READ_PARAMETER(reader, double, target_area_min);
        READ_PARAMETER(reader, double, target_area_max);
        READ_PARAMETER(reader, double, shape_index);
        READ_PARAMETER(reader, double, neightbor_alignment_intensity);
        READ_PARAMETER(reader, double, shape_alignment_intensity);
        READ_PARAMETER(reader, double, protrusion_alignment_intensity);
        READ_PARAMETER(reader, double, random_polarity_delta);
        READ_PARAMETER(reader, unsigned, correlation_sample_bucket_count);

        RandomNumberGenerator::Instance()->Reseed(random_seed);

        TruncatedNormal truncated_normal{ target_area_parent_normal_mean,
                                          target_area_parent_normal_variance,
                                          target_area_min,
                                          target_area_max };
        truncated_normal.Generate(cell_count);

        MyVoronoiVertexMeshGenerator generator{ truncated_normal.GetGenerated() };
        MutableVertexMesh<2, 2>* mesh = generator.GetMesh();
        if (truncated_normal.GetGenerated().size() != mesh->GetNumElements())
        {
            EXCEPTION("incompatible");
        }

        std::vector<CellPtr> cells{};
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<2> cell_population{ *mesh, cells };

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
        for (auto iter = mesh->GetNodeIteratorBegin(); iter != mesh->GetNodeIteratorEnd(); ++iter)
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
        p_nagai_honda_force->SetCellExpansionTerm(nagai_honda_cell_expansion_parameter);
        simulator.AddForce(p_nagai_honda_force);

        std::string const align_polarity{ "align polarity theta" };
        std::string const random_polarity{ "random polarity theta" };

        boost::shared_ptr<PolarityForce> align_polarity_force = boost::make_shared<PolarityForce>(align_propelling_parameter, align_polarity);
        simulator.AddForce(align_polarity_force);

        boost::shared_ptr<PolarityForce> random_polarity_force = boost::make_shared<PolarityForce>(random_propelling_parameter, random_polarity);
        simulator.AddForce(random_polarity_force);

        boost::shared_ptr<FixedTargetAreaModifier> fixed_target_area_modifier{};
        if (uniform_taret_area)
        {
            fixed_target_area_modifier = boost::make_shared<FixedTargetAreaModifier>(target_area_parent_normal_mean, shape_index);
            simulator.AddSimulationModifier(fixed_target_area_modifier);
        }
        else
        {
            fixed_target_area_modifier = boost::make_shared<FixedTargetAreaModifier>(truncated_normal, shape_index);
            simulator.AddSimulationModifier(fixed_target_area_modifier);
        }

        boost::shared_ptr<PolarityModifier> polarity_modifier = boost::make_shared<PolarityModifier>(neightbor_alignment_intensity,
                                                                                                     shape_alignment_intensity,
                                                                                                     protrusion_alignment_intensity,
                                                                                                     simulator.GetDt(),
                                                                                                     align_polarity);
        simulator.AddSimulationModifier(polarity_modifier);

        boost::shared_ptr<SimplePolarityModifier> simple_polarity_modifier = boost::make_shared<SimplePolarityModifier>(random_polarity_delta,
                                                                                                                        random_polarity);
        simulator.AddSimulationModifier(simple_polarity_modifier);

        auto p_vertex_velocity_correlation_tracker = boost::make_shared<VertexVelocityCorrelationTrackingModifier>(box_diameter, correlation_sample_bucket_count, file_finder.GetLeafName() + "VertexVelocityCorrelation");
        simulator.AddSimulationModifier(p_vertex_velocity_correlation_tracker);

        auto p_vertex_displacement_tracker = boost::make_shared<VertexDisplacementTrackingModifier>(simulator.GetDt(), file_finder.GetLeafName() + "VertexDisplacement");
        simulator.AddSimulationModifier(p_vertex_displacement_tracker);

        simulator.Solve();

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
    catch (Exception& exp)
    {
        TS_FAIL(exp.GetMessage());
    }
}

class TestFinal : public AbstractCellBasedTestSuite
{
public:
    // void TestNonuniform_30_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_00"); }
    // void TestNonuniform_30_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_01"); }
    // void TestNonuniform_30_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_02"); }
    // void TestNonuniform_30_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_03"); }
    // void TestNonuniform_30_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_04"); }
    // void TestNonuniform_30_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_05"); }
    // void TestNonuniform_30_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_06"); }
    // void TestNonuniform_30_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_07"); }
    // void TestNonuniform_30_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_08"); }
    // void TestNonuniform_30_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_09"); }
    // void TestNonuniform_30_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_30_10"); }

    // void TestNonuniform_31_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_00"); }
    // void TestNonuniform_31_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_01"); }
    // void TestNonuniform_31_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_02"); }
    // void TestNonuniform_31_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_03"); }
    // void TestNonuniform_31_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_04"); }
    // void TestNonuniform_31_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_05"); }
    // void TestNonuniform_31_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_06"); }
    // void TestNonuniform_31_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_07"); }
    // void TestNonuniform_31_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_08"); }
    // void TestNonuniform_31_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_09"); }
    // void TestNonuniform_31_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_31_10"); }

    // void TestNonuniform_32_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_00"); }
    // void TestNonuniform_32_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_01"); }
    // void TestNonuniform_32_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_02"); }
    // void TestNonuniform_32_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_03"); }
    // void TestNonuniform_32_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_04"); }
    // void TestNonuniform_32_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_05"); }
    // void TestNonuniform_32_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_06"); }
    // void TestNonuniform_32_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_07"); }
    // void TestNonuniform_32_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_08"); }
    // void TestNonuniform_32_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_09"); }
    // void TestNonuniform_32_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_32_10"); }

    // void TestNonuniform_33_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_00"); }
    // void TestNonuniform_33_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_01"); }
    // void TestNonuniform_33_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_02"); }
    // void TestNonuniform_33_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_03"); }
    // void TestNonuniform_33_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_04"); }
    // void TestNonuniform_33_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_05"); }
    // void TestNonuniform_33_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_06"); }
    // void TestNonuniform_33_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_07"); }
    // void TestNonuniform_33_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_08"); }
    // void TestNonuniform_33_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_09"); }
    // void TestNonuniform_33_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_33_10"); }

    // void TestNonuniform_34_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_00"); }
    // void TestNonuniform_34_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_01"); }
    // void TestNonuniform_34_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_02"); }
    // void TestNonuniform_34_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_03"); }
    // void TestNonuniform_34_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_04"); }
    // void TestNonuniform_34_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_05"); }
    // void TestNonuniform_34_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_06"); }
    // void TestNonuniform_34_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_07"); }
    // void TestNonuniform_34_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_08"); }
    // void TestNonuniform_34_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_09"); }
    // void TestNonuniform_34_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_34_10"); }

    // void TestNonuniform_35_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_00"); }
    // void TestNonuniform_35_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_01"); }
    // void TestNonuniform_35_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_02"); }
    // void TestNonuniform_35_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_03"); }
    // void TestNonuniform_35_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_04"); }
    // void TestNonuniform_35_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_05"); }
    // void TestNonuniform_35_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_06"); }
    // void TestNonuniform_35_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_07"); }
    // void TestNonuniform_35_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_08"); }
    // void TestNonuniform_35_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_09"); }
    // void TestNonuniform_35_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_35_10"); }

    // void TestNonuniform_36_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_00"); }
    // void TestNonuniform_36_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_01"); }
    // void TestNonuniform_36_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_02"); }
    // void TestNonuniform_36_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_03"); }
    // void TestNonuniform_36_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_04"); }
    // void TestNonuniform_36_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_05"); }
    // void TestNonuniform_36_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_06"); }
    // void TestNonuniform_36_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_07"); }
    // void TestNonuniform_36_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_08"); }
    // void TestNonuniform_36_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_09"); }
    // void TestNonuniform_36_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_36_10"); }

    // void TestNonuniform_37_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_00"); }
    // void TestNonuniform_37_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_01"); }
    // void TestNonuniform_37_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_02"); }
    // void TestNonuniform_37_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_03"); }
    // void TestNonuniform_37_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_04"); }
    // void TestNonuniform_37_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_05"); }
    // void TestNonuniform_37_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_06"); }
    // void TestNonuniform_37_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_07"); }
    // void TestNonuniform_37_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_08"); }
    // void TestNonuniform_37_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_09"); }
    // void TestNonuniform_37_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_37_10"); }

    // void TestNonuniform_38_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_00"); }
    // void TestNonuniform_38_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_01"); }
    // void TestNonuniform_38_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_02"); }
    // void TestNonuniform_38_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_03"); }
    // void TestNonuniform_38_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_04"); }
    // void TestNonuniform_38_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_05"); }
    // void TestNonuniform_38_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_06"); }
    // void TestNonuniform_38_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_07"); }
    // void TestNonuniform_38_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_08"); }
    // void TestNonuniform_38_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_09"); }
    // void TestNonuniform_38_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_38_10"); }

    // void TestNonuniform_39_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_00"); }
    // void TestNonuniform_39_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_01"); }
    // void TestNonuniform_39_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_02"); }
    // void TestNonuniform_39_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_03"); }
    // void TestNonuniform_39_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_04"); }
    // void TestNonuniform_39_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_05"); }
    // void TestNonuniform_39_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_06"); }
    // void TestNonuniform_39_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_07"); }
    // void TestNonuniform_39_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_08"); }
    // void TestNonuniform_39_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_09"); }
    // void TestNonuniform_39_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_39_10"); }

    // void TestNonuniform_40_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_00"); }
    // void TestNonuniform_40_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_01"); }
    // void TestNonuniform_40_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_02"); }
    // void TestNonuniform_40_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_03"); }
    // void TestNonuniform_40_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_04"); }
    // void TestNonuniform_40_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_05"); }
    // void TestNonuniform_40_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_06"); }
    // void TestNonuniform_40_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_07"); }
    void TestNonuniform_40_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_08"); }
    // void TestNonuniform_40_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_09"); }
    // void TestNonuniform_40_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalNonuniform_40_10"); }

public:
    // void TestUniform_30_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_00"); }
    // void TestUniform_30_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_01"); }
    // void TestUniform_30_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_02"); }
    // void TestUniform_30_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_03"); }
    // void TestUniform_30_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_04"); }
    // void TestUniform_30_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_05"); }
    // void TestUniform_30_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_06"); }
    // void TestUniform_30_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_07"); }
    // void TestUniform_30_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_08"); }
    // void TestUniform_30_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_09"); }
    // void TestUniform_30_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_30_10"); }

    // void TestUniform_31_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_00"); }
    // void TestUniform_31_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_01"); }
    // void TestUniform_31_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_02"); }
    // void TestUniform_31_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_03"); }
    // void TestUniform_31_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_04"); }
    // void TestUniform_31_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_05"); }
    // void TestUniform_31_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_06"); }
    // void TestUniform_31_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_07"); }
    // void TestUniform_31_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_08"); }
    // void TestUniform_31_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_09"); }
    // void TestUniform_31_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_31_10"); }

    // void TestUniform_32_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_00"); }
    // void TestUniform_32_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_01"); }
    // void TestUniform_32_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_02"); }
    // void TestUniform_32_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_03"); }
    // void TestUniform_32_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_04"); }
    // void TestUniform_32_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_05"); }
    // void TestUniform_32_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_06"); }
    // void TestUniform_32_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_07"); }
    // void TestUniform_32_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_08"); }
    // void TestUniform_32_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_09"); }
    // void TestUniform_32_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_32_10"); }

    // void TestUniform_33_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_00"); }
    // void TestUniform_33_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_01"); }
    // void TestUniform_33_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_02"); }
    // void TestUniform_33_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_03"); }
    // void TestUniform_33_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_04"); }
    // void TestUniform_33_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_05"); }
    // void TestUniform_33_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_06"); }
    // void TestUniform_33_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_07"); }
    // void TestUniform_33_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_08"); }
    // void TestUniform_33_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_09"); }
    // void TestUniform_33_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_33_10"); }

    // void TestUniform_34_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_00"); }
    // void TestUniform_34_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_01"); }
    // void TestUniform_34_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_02"); }
    // void TestUniform_34_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_03"); }
    // void TestUniform_34_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_04"); }
    // void TestUniform_34_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_05"); }
    // void TestUniform_34_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_06"); }
    // void TestUniform_34_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_07"); }
    // void TestUniform_34_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_08"); }
    // void TestUniform_34_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_09"); }
    // void TestUniform_34_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_34_10"); }

    // void TestUniform_35_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_00"); }
    // void TestUniform_35_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_01"); }
    // void TestUniform_35_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_02"); }
    // void TestUniform_35_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_03"); }
    // void TestUniform_35_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_04"); }
    // void TestUniform_35_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_05"); }
    // void TestUniform_35_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_06"); }
    // void TestUniform_35_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_07"); }
    // void TestUniform_35_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_08"); }
    // void TestUniform_35_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_09"); }
    // void TestUniform_35_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_35_10"); }

    // void TestUniform_36_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_00"); }
    // void TestUniform_36_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_01"); }
    // void TestUniform_36_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_02"); }
    // void TestUniform_36_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_03"); }
    // void TestUniform_36_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_04"); }
    // void TestUniform_36_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_05"); }
    // void TestUniform_36_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_06"); }
    // void TestUniform_36_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_07"); }
    // void TestUniform_36_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_08"); }
    // void TestUniform_36_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_09"); }
    // void TestUniform_36_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_36_10"); }

    // void TestUniform_37_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_00"); }
    // void TestUniform_37_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_01"); }
    // void TestUniform_37_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_02"); }
    // void TestUniform_37_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_03"); }
    // void TestUniform_37_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_04"); }
    // void TestUniform_37_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_05"); }
    // void TestUniform_37_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_06"); }
    // void TestUniform_37_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_07"); }
    // void TestUniform_37_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_08"); }
    // void TestUniform_37_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_09"); }
    // void TestUniform_37_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_37_10"); }

    // void TestUniform_38_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_00"); }
    // void TestUniform_38_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_01"); }
    // void TestUniform_38_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_02"); }
    // void TestUniform_38_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_03"); }
    // void TestUniform_38_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_04"); }
    // void TestUniform_38_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_05"); }
    // void TestUniform_38_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_06"); }
    // void TestUniform_38_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_07"); }
    // void TestUniform_38_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_08"); }
    // void TestUniform_38_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_09"); }
    // void TestUniform_38_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_38_10"); }

    // void TestUniform_39_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_00"); }
    // void TestUniform_39_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_01"); }
    // void TestUniform_39_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_02"); }
    // void TestUniform_39_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_03"); }
    // void TestUniform_39_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_04"); }
    // void TestUniform_39_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_05"); }
    // void TestUniform_39_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_06"); }
    // void TestUniform_39_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_07"); }
    // void TestUniform_39_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_08"); }
    // void TestUniform_39_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_09"); }
    // void TestUniform_39_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_39_10"); }

    // void TestUniform_40_00() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_00"); }
    // void TestUniform_40_01() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_01"); }
    // void TestUniform_40_02() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_02"); }
    // void TestUniform_40_03() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_03"); }
    // void TestUniform_40_04() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_04"); }
    // void TestUniform_40_05() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_05"); }
    // void TestUniform_40_06() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_06"); }
    // void TestUniform_40_07() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_07"); }
    // void TestUniform_40_08() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_08"); }
    // void TestUniform_40_09() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_09"); }
    // void TestUniform_40_10() { ParameterizedTest("cell_based/test/pozero/data/FinalBase", "cell_based/test/pozero/data/FinalUniform_40_10"); }
};

#endif
