#ifndef SIMPLEPOLARITYMODIFIER_HPP_
#define SIMPLEPOLARITYMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

class SimplePolarityModifier : public AbstractCellBasedSimulationModifier<2>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialization(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
        archive & mPolarityDelta;
    }

    double mPolarityDelta;

    std::string mPolarityName;

public:
    SimplePolarityModifier(double polarity_delta, std::string polarity_name);

    virtual ~SimplePolarityModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double GetPolarityDelta() const;

    void SetPolarityDelta(double newParam);
};

#endif /*SIMPLEPOLARITYMODIFIER_HPP_*/
