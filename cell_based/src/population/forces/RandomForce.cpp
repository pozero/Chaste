#include "RandomForce.hpp"

RandomForce::RandomForce(double scale) : mScale(scale)
{
}

RandomForce::~RandomForce()
{
}

void RandomForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("RandomForce expects to be used with a VertexBasedCellPopulation");
    }
    RandomNumberGenerator* rng = RandomNumberGenerator::Instance();
    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned const num_nodes = p_cell_population->GetNumNodes();
    for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        c_vector<double, 2> contribution = zero_vector<double>(2);
        contribution[0] = rng->ranf();
        contribution[1] = rng->ranf();
        contribution = contribution / norm_2(contribution);
        contribution = mScale * contribution;
        p_cell_population->GetNode(node_idx)->AddAppliedForceContribution(contribution);
    }
}

void RandomForce::OutputForceParameters(out_stream& rParamsFile)
{
}
