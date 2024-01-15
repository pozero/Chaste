# TruncatedNormalDistribution

We need to construct a non-zero normal-like distribution to generate cells' target areas. 

## Sampling

The sampling is intuitive, just sample from the parent normal distribution then reject all value outside of the scope, specifically, less than or equal zero.

```c++
RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
double target_area = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
while (target_area <= 0.0) {
    target_area = p_gen->NormalRandomDeviate(mParentNormalMean, mParentNormalVariance);
}
```

## Quantile

Calculation of quantile is a little bit more complicated. We know the CDF of truncated normal distribution is

$$
\Psi(\overline{\mu}, \overline{\sigma}, a, b; x) = \frac{\Phi(\overline{\mu}, \overline{\sigma}; x) - \Phi(\overline{\mu}, \overline{\sigma}; a)}{ \Phi(\overline{\mu}, \overline{\sigma}; b) - \Phi(\overline{\mu}, \overline{\sigma}; a)}
$$

And the CDF of normal distribution is

$$
\Phi(\overline{\mu}, \overline{\sigma}; x) = \frac12 \left[ 1 + \text{erf}\left( \frac{x - \mu}{\sqrt{2}\sigma} \right) \right]
$$

Let the desired quantile be $p$, we get the equation

$$
\begin {align}
p &= \frac{\Phi(\overline{\mu}, \overline{\sigma}; x) - \Phi_a}{\Phi_b - \Phi_a}\\
p &= \frac{\frac12 \left[1 + \text{erf}\left( \frac{x - \mu}{\sqrt{2}\sigma} \right) \right] - \Phi_a}{\Phi_b - \Phi_a}\\
2p(\Phi_b - \Phi_a) + 2\Phi_a - 1 &= \text{erf}\left( \frac{x - \mu}{\sqrt{2}\sigma} \right)\\
\text{erf}^{-1}(2p(\Phi_b - \Phi_a) + 2\Phi_a - 1) &= \frac{x - \mu}{\sqrt{2}\sigma}\\
\sqrt{2}\sigma\text{erf}^{-1}(2p(\Phi_b - \Phi_a) + 2\Phi_a - 1) + \mu &= x\\
\end {align}
$$

