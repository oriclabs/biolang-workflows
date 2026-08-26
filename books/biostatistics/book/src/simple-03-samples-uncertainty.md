# Start Simple 3: Samples, Standard Error, and Confidence Intervals

> **In one sentence:** SD describes differences among observations; standard
> error describes how much an estimate would vary across repeated samples.

## Population, Sample, Parameter, Statistic

Suppose we want the mean systolic blood pressure of every eligible adult in a
region. Measuring everyone is impractical, so we measure a sample.

```text
Target population                    Random sample
o o o o o o o o o o                 o   o  o
 o o o o o o o o o       ---->        o  o
o o o o o o o o o o                 o   o  o

unknown population mean             observed sample mean
      parameter                           statistic
```

A **parameter** describes the population. A **statistic** is calculated from a
sample and estimates that parameter. Another sample would usually give a
slightly different statistic. That ordinary movement is sampling variation.

## SD Is Not SE

| Quantity | Question | Usually changes when n grows? |
|---|---|---|
| SD | How different are individual observations? | Not simply because n grows |
| SE of the mean | How variable would sample means be? | Usually decreases as `1 / sqrt(n)` |
| Confidence interval | Which parameter values remain reasonably compatible with the data and model? | Usually narrows as information grows |

For independent observations under the usual sample-mean model:

```text
SE = sample SD / sqrt(sample size)
```

Four times as many independent observations gives roughly half the SE, not one
quarter. Repeated measurements from the same person are not automatically four
independent observations.

## Calculate a Mean and Interval

```biolang
let pressure_values = [121.0, 128.0, 119.0, 133.0, 126.0,
                       124.0, 130.0, 117.0, 129.0, 125.0]

let n = len(pressure_values)
let estimate = mean(pressure_values)
let se = stdev(pressure_values) / sqrt(n)

# 2.262 is the two-sided 95% t critical value for 9 degrees of freedom.
let interval = [estimate - 2.262 * se, estimate + 2.262 * se]

{
  sample_mean: estimate,
  sample_sd: stdev(pressure_values),
  standard_error: se,
  confidence_interval_95: interval
}
```

The multiplier is about 2 here but is not always exactly 2. It depends on the
model, confidence level, sample size, and degrees of freedom.

## What a 95% Confidence Interval Means

A 95% confidence procedure is designed so that, over many compatible repeated
samples, about 95% of the intervals it constructs cover the fixed population
parameter.

For one observed interval, a useful reporting interpretation is:

> Under the sampling and analysis assumptions, values inside this interval are
> reasonably compatible with the observed data; values outside are less
> compatible.

It does not mean that 95% of individual observations lie inside the interval.
It also does not repair selection bias or poor measurement.

## Sampling Methods Matter

- **Simple random:** each eligible unit has a known equal selection chance.
- **Stratified:** sample separately within important groups.
- **Systematic:** choose a random start and then every kth unit, after checking
  that list order has no matching pattern.
- **Cluster:** sample groups such as clinics, then observe units within them.
- **Convenience:** take what is easiest; useful for exploration but vulnerable
  to poor generalization.

Clustered and stratified designs may require design-aware SE calculations. A
large convenience sample can be precisely wrong about the target population.

## A Reporting Sentence

> In 10 sampled adults, mean systolic pressure was reported with its SD to
> describe individual variation and a 95% confidence interval to describe
> uncertainty in the population mean.

> **Do not conclude:** that a narrow interval proves the sample is
> representative. Precision and bias are different problems.

## Quick Check

1. Which quantity describes individual variation: SD or SE?
2. If independent sample size grows from 25 to 100, approximately what happens
   to the SE?
3. Can 1,000 measurements from one dish replace 1,000 independent dishes?

<details>
<summary>Answers</summary>

1. SD.
2. It is approximately halved, if the underlying variability and design stay
   comparable.
3. No. Those measurements share one experimental unit and are dependent.

</details>

**Next:** [Probability, Tests, and p-Values](simple-04-probability-tests.md).
For more depth, see [Day 5](day-05.md) and [Day 6](day-06.md).
