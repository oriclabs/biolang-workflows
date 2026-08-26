# Start Simple 2: Spread, Shape, Outliers, and Logs

> **In one sentence:** Spread describes how far observations differ; shape
> shows where those differences occur.

## Same Centre, Different Data

Both groups below have mean 10.

```text
Tight:    8  9  10  11  12        small spread
Wide:     0  5  10  15  20        large spread
                   ^
                same mean
```

The mean alone cannot tell these groups apart. Spread supplies the missing
information.

## Measures of Spread

| Spread | What it describes | Pair it with |
|---|---|---|
| **Range** | Lowest to highest; very easy but controlled by two values | A quick first look |
| **IQR** | Width of the middle 50% | Median |
| **MAD** | Typical absolute distance from the median | Median |
| **Variance** | Average squared distance from the mean, with sample correction | Models and calculations |
| **SD** | Square root of variance, returned to the original measurement unit | Mean |

Variance is not the "whole spread" while SD is the left and right half. They
describe the same variability on different scales. If SD is 3 mg/L, variance
is 9 `(mg/L)^2`. SD is easier to interpret because its unit matches the data.

```biolang
let tight_values = [8.0, 9.0, 10.0, 11.0, 12.0]
let wide_values = [0.0, 5.0, 10.0, 15.0, 20.0]

let spread_comparison = {
  tight: {
    mean: mean(tight_values),
    sd: stdev(tight_values),
    variance: variance(tight_values),
    iqr: quantile(tight_values, 0.75) - quantile(tight_values, 0.25)
  },
  wide: {
    mean: mean(wide_values),
    sd: stdev(wide_values),
    variance: variance(wide_values),
    iqr: quantile(wide_values, 0.75) - quantile(wide_values, 0.25)
  }
}

spread_comparison
```

## What Does Mean Plus or Minus SD Mean?

For an approximately bell-shaped distribution:

```text
                 about 68%
             <-------------->
        about 95% within mean +/- 2 SD
   about 99.7% within mean +/- 3 SD

----|---------|---------|---------|---------|---------|----
  -3 SD     -2 SD     -1 SD      mean     +1 SD     +2 SD
```

These percentages are a property of a normal model, not a rule for every
dataset. Strong skew, several peaks, bounds, and many zeros can make the rule
misleading.

An observation beyond 3 SD is unusual under a normal model. It is not
automatically an error and should not be deleted automatically.

## Read Shape Before Choosing a Summary

```text
Symmetric:       ▁▃▆█▆▃▁       mean and median are close
Right-skewed:    █▆▄▃▂▁        mean is pulled toward the long right tail
Two groups:      ▂▆█▃  ▂▇█▃    one centre may describe neither group
Many zeros:      █          ▂▁  absence and positive values may need separate thought
```

A histogram is a visual clue, not an automatic test selector.

```biolang
let count_values = [0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 8.0, 20.0]
histogram(count_values, {bins: 8, title: "Right-skewed counts"})
```

## When Does a Log Help?

A log transformation is useful when values are positive and ratios or fold
changes are scientifically meaningful. It can compress a long right tail:

```text
Raw:     1 ---- 10 ------------------------------ 100
log10:   0 ----  1 ------------------------------- 2

Equal distance after logging means equal multiplication before logging.
```

Use `log2(x + 1)` for non-negative counts only when the added 1 and log scale
match the analysis goal. Do not log negative values, silently add arbitrary
constants, or log data merely to obtain a desired p-value.

## A Reporting Sentence

> Both groups had mean 10, but the wide group had much greater SD and IQR;
> therefore the mean alone did not describe their difference.

> **Do not conclude:** that SD measures uncertainty in the mean. SD describes
> observations. Standard error, introduced next, describes uncertainty in an
> estimate.

## Quick Check

1. Which spread usually accompanies the median?
2. If SD is 4 mm, what is the variance and what unit does it use?
3. What biological meaning should justify a log transformation?

<details>
<summary>Answers</summary>

1. IQR or MAD.
2. Variance is 16 mm-squared.
3. Multiplicative differences, ratios, or fold changes should be meaningful;
   the data must also be valid on the chosen log scale.

</details>

**Next:** [Samples, Standard Error, and Confidence Intervals](simple-03-samples-uncertainty.md).
