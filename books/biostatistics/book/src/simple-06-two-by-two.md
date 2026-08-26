# Start Simple 6: One 2x2 Table, Many Useful Measures

> **In one sentence:** A 2x2 table can describe absolute risk, risk difference,
> risk ratio, odds ratio, and number needed to treat - but the study design
> determines which measures are valid.

## Build the Table First

Suppose a cohort records an exposure and later observes disease:

| | Disease | No disease | Total |
|---|---:|---:|---:|
| Exposed | `a = 30` | `b = 70` | 100 |
| Unexposed | `c = 10` | `d = 90` | 100 |

```text
Risk among exposed   = a / (a + b) = 30 / 100
Risk among unexposed = c / (c + d) = 10 / 100
```

The two risks are the foundation. Differences ask an additive question;
ratios ask a multiplicative question.

## Calculate the Main Measures

```biolang
let a = 30.0
let b = 70.0
let c = 10.0
let d = 90.0

let risk_exposed = a / (a + b)
let risk_unexposed = c / (c + d)
let risk_difference = risk_exposed - risk_unexposed
let risk_ratio = risk_exposed / risk_unexposed
let odds_ratio = (a * d) / (b * c)
let number_needed_to_harm = 1.0 / abs(risk_difference)
let exact_test = fisher_exact(a, b, c, d)

{
  risk_exposed: risk_exposed,
  risk_unexposed: risk_unexposed,
  risk_difference: risk_difference,
  risk_ratio: risk_ratio,
  odds_ratio: odds_ratio,
  number_needed_to_harm: number_needed_to_harm,
  fisher_p_value: exact_test.p_value
}
```

For these data, exposed risk is 30%, unexposed risk is 10%, the risk
difference is 20 percentage points, and the risk ratio is 3. These statements
describe different aspects of the same table.

## What Each Measure Says

| Measure | Question | Null value |
|---|---|---:|
| Risk difference | How many additional cases occur per person exposed? | 0 |
| Risk ratio | How many times as large is the exposed risk? | 1 |
| Odds ratio | How many times as large are the exposed odds? | 1 |
| NNT or NNH | How many people need treatment or exposure for one additional outcome over the stated time? | No simple null |

When treatment lowers risk, use the absolute risk reduction to calculate NNT.
When exposure raises harmful risk, the same reciprocal is often called NNH.
Always state the follow-up period.

## Odds Are Not Risk

If risk is `p`, odds are `p / (1 - p)`.

```text
Risk 10%  -> odds 1 to 9
Risk 50%  -> odds 1 to 1
Risk 80%  -> odds 4 to 1
```

When outcomes are common, an odds ratio can look much farther from 1 than the
risk ratio. Do not translate it casually into "times more likely."

## Study Design Controls Interpretation

- In a cohort or suitable randomized trial, risks and risk ratios can usually
  be calculated from follow-up counts.
- In a conventional case-control study, the investigator chooses the number of
  cases and controls. Those sampled counts do not estimate population risk.
  The odds ratio remains the standard association measure.
- In any observational design, association may reflect confounding or bias.

The p-value from `fisher_exact()` concerns evidence against no association. It
does not replace the effect measures or establish causality.

## A Reporting Sentence

> Over the stated follow-up, disease risk was 30% in exposed participants and
> 10% in unexposed participants: a risk difference of 20 percentage points and
> a risk ratio of 3.0.

> **Do not conclude:** "three times the risk" without also reporting the
> absolute risks. A large ratio can describe a very small absolute change.

## Quick Check

1. What is the null value for a risk ratio?
2. Which measure directly gives additional cases per exposed person?
3. Why can a case-control table usually not provide population risk?

<details>
<summary>Answers</summary>

1. One.
2. Risk difference.
3. The researcher fixes or samples the numbers of cases and controls; their
   proportions are not the population's natural outcome proportions.

</details>

**Next:** [Sampling and Study Designs](simple-07-study-designs.md). For more
categorical analysis, see [Day 11](day-11.md) and [Day 19](day-19.md).
