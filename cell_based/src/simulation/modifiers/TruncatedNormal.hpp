#ifndef TRUNCATEDNORMAL_HPP_
#define TRUNCATEDNORMAL_HPP_
#include "RandomNumberGenerator.hpp"

class TruncatedNormal
{
private:
    double mParentNormalMean;

    double mParentNormalVariance;

    double mMin;

    double mMax;

    std::vector<double> mGenerated{};

public:
    TruncatedNormal(double parent_mean, double parent_variance, double min, double max)
            : mParentNormalMean(parent_mean), mParentNormalVariance(parent_variance), mMin(min), mMax(max) {}

    double Generate() const
    {

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
        while (val < mMin || val > mMax)
        {
            val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
        }
        return val;
    }

    void Generate(uint32_t count)
    {
        for (uint32_t i = 0; i < count; ++i)
        {
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            double val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
            while (val < mMin || val > mMax)
            {
                val = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
            }
            mGenerated.push_back(val);
        }
    }

    std::vector<double> const& GetGenerated() const
    {
        return mGenerated;
    }

    std::vector<double>& GetGenerated()
    {
        return mGenerated;
    }

    double ParentNormalCDF(double val) const
    {
        return 0.5 * (1.0 + boost::math::erf((val - mParentNormalMean) / (sqrt(2.0) * mParentNormalVariance)));
    }

    double TruncatedNormalQuantile(double phi_a, double phi_b, double percentage) const
    {
        return sqrt(2.0) * mParentNormalVariance * boost::math::erf_inv(2.0 * percentage * (phi_b - phi_a) + 2.0 * phi_a - 1.0) + mParentNormalMean;
    }

    double GetMean() const
    {
        return mParentNormalMean;
    }

    double GetVariance() const
    {
        return mParentNormalVariance;
    }

    double GetMin() const
    {
        return mMin;
    }

    double GetMax() const
    {
        return mMax;
    }
};
#endif
