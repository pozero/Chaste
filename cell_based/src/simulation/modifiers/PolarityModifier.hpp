#ifndef POLARITYMODIFIER_HPP_
#define POLARITYMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

class PolarityModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialization(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
        archive & mNeighborAlignmentIntensity;
    }

    double mNeighborAlignmentIntensity;

    double mShapeAlignmentIntensity;

    // double mProtrusionAlignmentIntensity;

    double mDt;

    // std::vector<c_vector<double, 2> > mCellOutwradVectorSum{};

public:
    PolarityModifier();

    virtual ~PolarityModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double GetNeighborAlignmentIntensity() const;

    double GetShapeAlignmentIntensity() const;

    // double GetProtrusionAlignmentIntensity() const;

    void SetNeighborAlignmentIntensity(double newParam);

    void SetShapeAlignmentIntensity(double newParam);

    // void SetProtrusionAlignmentIntensity(double newParam);

    void SetDt(double dt);
};

#endif
