#ifndef MYHONGDAFORCE_HPP_
#define MYHONGDAFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

template <unsigned DIM>
class MyNagaiHondaForce : public AbstractForce<DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mNagaiHondaDeformationEnergyParameter;
        archive & mNagaiHondaMembraneSurfaceEnergyParameter;
    }

    double mNagaiHondaDeformationEnergyParameter;

    double mNagaiHondaMembraneSurfaceEnergyParameter;

    double mCellExpansionTerm;

public:
    MyNagaiHondaForce();

    virtual ~MyNagaiHondaForce();

    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    double GetDeformationEnergyParameter() const;

    double GetMembraneSurfaceEnergyParameter() const;

    void SetDeformationEnergyParameter(double newParam);

    void SetMembraneSurfaceEnergyParameter(double newParam);

    void SetCellExpansionTerm(double newParam);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNagaiHondaForce)

#endif /*MYHONGDAFORCE_HPP_*/
