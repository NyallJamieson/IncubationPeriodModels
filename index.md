# IncubationPeriodModels

`IncubationPeriodModels` implements flexible Burr-family distributions
for modelling incubation periods and other epidemiological delay
distributions.

The package provides density, distribution, quantile, and random
generation functions for several Burr-family distributions and a derived
Burr-like distribution. These models offer flexible parametric
representations of positive-valued event times and can accommodate a
wide range of skewness, tail behaviour, and hazard shapes.

Such flexibility makes Burr distributions useful for modelling
incubation periods and other delay distributions encountered in
epidemiology where heavy tails and asymmetric shapes are common.

## Implemented distributions

The package currently provides implementations of:

- Burr Type III
- Burr Type X
- Burr Type XII
- Derived Burr distribution

All distributions follow the standard R distribution interface:

| Function | Purpose                 |
|----------|-------------------------|
| `d*`     | density                 |
| `p*`     | cumulative distribution |
| `q*`     | quantile                |
| `r*`     | random generation       |

## Installation

Install the development version from GitHub:

\`\`\`r \# install.packages(“pak”)
pak::pak(“NyallJamieson/IncubationPeriodModels”)
