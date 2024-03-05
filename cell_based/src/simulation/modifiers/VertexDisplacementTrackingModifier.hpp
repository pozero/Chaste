#ifndef VERTEXDISPLACEMENTTRACKINGMODIFIER_HPP_
#define VERTEXDISPLACEMENTTRACKINGMODIFIER_HPP_

#include <vector>

#include "AbstractCellBasedSimulationModifier.hpp"

class VertexDisplacementTrackingModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mDt;

    std::string mOutputFileName;

    std::vector<double> mMSD;

    std::ofstream mOutputFile;

public:
    VertexDisplacementTrackingModifier() = delete;

    VertexDisplacementTrackingModifier(double dt, std::string const& output_file);

    virtual ~VertexDisplacementTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#endif
