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

    std::string mPolarityName;

public:
    PolarityForce(double self_propelling, std::string polarity_name);

    virtual ~PolarityForce();

    virtual void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#endif /*POLARITYFORCE_HPP_*/
