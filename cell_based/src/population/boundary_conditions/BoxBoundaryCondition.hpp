#ifndef BOXBOUNDARYCONDITION_HPP_
#define BOXBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

class BoxBoundaryCondition : public AbstractCellPopulationBoundaryCondition<2>
{
private:
    double m_top;
    double m_bottom;
    double m_left;
    double m_right;

public:
    BoxBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation, double top, double bottom, double left, double right);

    virtual void ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations);

    virtual bool VerifyBoundaryCondition();

    virtual void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#endif
