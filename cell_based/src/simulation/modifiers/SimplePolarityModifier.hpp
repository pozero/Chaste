#ifndef SIMPLEPOLARITYMODIFIER_HPP_
#define SIMPLEPOLARITYMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

template <unsigned DIM>
class SimplePolarityModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialization(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
    }

    double mPolarityDelta;

public:
    SimplePolarityModifier();

    virtual ~SimplePolarityModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double GetPolarityDelta() const;

    void SetPolarityDelta(double newParam);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimplePolarityModifier)

#endif /*SIMPLEPOLARITYMODIFIER_HPP_*/
