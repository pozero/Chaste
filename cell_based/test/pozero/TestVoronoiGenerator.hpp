#ifndef TESTVORONOIGENERATOR_HPP_
#define TESTVORONOIGENERATOR_HPP_

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

class TestVoronoiGenerator : public AbstractCellBasedTestSuite
{
public:
    // void TestBase() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorBase"); }

    void TestNonUniform1() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform1"); }
    // void TestNonUniform2() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform2"); }
    // void TestNonUniform3() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform3"); }
    // void TestNonUniform4() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform4"); }
    // void TestNonUniform5() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform5"); }
    // void TestNonUniform6() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform6"); }
    // void TestNonUniform7() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform7"); }
    // void TestNonUniform8() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform8"); }
    // void TestNonUniform9() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform9"); }
    // void TestNonUniform10() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorNonUniform10"); }

    void TestUniform1() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform1"); }
    // void TestUniform2() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform2"); }
    // void TestUniform3() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform3"); }
    // void TestUniform4() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform4"); }
    // void TestUniform5() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform5"); }
    // void TestUniform6() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform6"); }
    // void TestUniform7() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform7"); }
    // void TestUniform8() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform8"); }
    // void TestUniform9() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform9"); }
    // void TestUniform10() { ParameterizedTest("cell_based/test/pozero/data/VoronoiGeneratorBase", "cell_based/test/pozero/data/VoronoiGeneratorUniform10"); }
};

#endif
