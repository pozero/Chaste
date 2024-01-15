#ifndef FIXEDTARGETAREAMODIFIER_HPP_
#define FIXEDTARGETAREAMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

#define FIXED_TARGET_AREA_SIZE_TAG "fixed target area size tag"

#define FIXED_TARGET_AREA_SMALL 1.0
#define FIXED_TARGET_AREA_MEDIUM 2.0
#define FIXED_TARGET_AREA_LARGE 3.0

template <unsigned DIM>
class FixedTargetAreaModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialization(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
        archive & mParentNormalMean;
        archive & mParentNormalVariance;
    }

    double mParentNormalMean;

    double mParentNormalVariance;

    double TruncatedNormal() const;

    double ParentNormalCDF(double val) const;

    double TruncatedNormalQuantile(double phi_a, double phi_b, double precentage) const;

public:
    FixedTargetAreaModifier();

    virtual ~FixedTargetAreaModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double GetParentNormalMean() const;

    double GetParentNormalVariance() const;

    void SetParentNormalMean(double newParam);

    void SetParentNormalVariance(double newParam);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedTargetAreaModifier)

#endif /*FIXEDTARGETAREAMODIFIER_HPP_*/
