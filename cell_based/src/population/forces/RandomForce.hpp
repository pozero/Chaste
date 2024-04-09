#ifndef RANDOMFORCE_HPP_
#define RANDOMFORCE_HPP_

#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "TruncatedNormal.hpp"

class RandomForce : public AbstractForce<2>
{
private:
    TruncatedNormal mScaleDist;

public:
    RandomForce(double scale_mean, double scale_variance, double scale_maximum);

    virtual ~RandomForce();

    virtual void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif
