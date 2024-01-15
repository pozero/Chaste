#ifndef POLARITYFORCE_HPP_
#define POLARITYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

class PolarityForce : public AbstractForce<2>
{

private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mSelfPropellingParameter;
    }

    double mSelfPropellingParameter;

public:
    PolarityForce();

    virtual ~PolarityForce();

    virtual void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    double GetSelfPropellingParameter() const;

    void SetSelfPropellingParameter(double newParam);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif /*POLARITYFORCE_HPP_*/
