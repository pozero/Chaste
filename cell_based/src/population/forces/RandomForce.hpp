#ifndef RANDOMFORCE_HPP_
#define RANDOMFORCE_HPP_

#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

class RandomForce : public AbstractForce<2>
{
private:
    double mScale;

public:
    RandomForce(double scale);

    virtual ~RandomForce();

    virtual void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif
