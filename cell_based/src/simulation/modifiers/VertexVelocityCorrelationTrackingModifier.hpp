#ifndef VERTEXVELOCITYCORRELATIONTRACKINGMODIFIER_HPP_
#define VERTEXVELOCITYCORRELATIONTRACKINGMODIFIER_HPP_

#include <map>
#include <fstream>

#include "AbstractCellBasedSimulationModifier.hpp"

class VertexVelocityCorrelationTrackingModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    double mBucketWidth;

    std::string mOutputFileName;

    std::map<double, double> mCorrelations;

    std::map<double, double> mDistances;

    std::map<double, int> mBucketCounts;

    std::ofstream mOutputFile;

public:
    VertexVelocityCorrelationTrackingModifier() = delete;

    VertexVelocityCorrelationTrackingModifier(double bucket_width, std::string const& output_file);

    virtual ~VertexVelocityCorrelationTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#endif
