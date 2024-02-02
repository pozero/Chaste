#ifndef INELASTICAL_CIRCLE_BOUNDARY_CONDITION_HPP_
#define INELASTICAL_CIRCLE_BOUNDARY_CONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

class InelasticalCircleBoundaryCondition : public AbstractCellPopulationBoundaryCondition<2>
{
private:
    c_vector<double, 2> m_circle_center;

    double m_circle_radius;

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
        archive & m_circle_center;
        archive & m_circle_radius;
    }

public:
    InelasticalCircleBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation,
                                       c_vector<double, 2> const& center,
                                       double radius);

    virtual void ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations);

    virtual bool VerifyBoundaryCondition();

    virtual void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#endif // INELASTICAL_CIRCLE_BOUNDARY_CONDITION_HPP_
