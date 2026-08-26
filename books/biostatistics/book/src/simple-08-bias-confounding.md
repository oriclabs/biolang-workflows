# Start Simple 8: Bias, Confounding, and Fair Comparisons

> **In one sentence:** Random error makes estimates wobble; bias moves them
> systematically; confounding mixes the effect of an exposure with another
> cause.

## Three Sources of Disagreement

```text
Random variation:   repeated estimates scatter around a target
Bias:               repeated estimates miss in one direction
Confounding:        groups differ in another cause of the outcome
```

Larger samples usually reduce random uncertainty. They do not automatically
remove bias or confounding. A million systematically mismeasured observations
can produce a very precise wrong answer.

## Common Forms of Bias

| Bias | What goes wrong | Example safeguard |
|---|---|---|
| Selection bias | Entry or retention relates to exposure and outcome | Define eligibility, recruitment, and follow-up clearly |
| Information bias | Exposure or outcome is measured differently | Use consistent, blinded measurement where possible |
| Recall bias | Memory differs between comparison groups | Prefer records or prospective measurement when suitable |
| Attrition bias | Loss to follow-up differs meaningfully | Report flow, reasons, and sensitivity analyses |
| Selective reporting | Only favorable outcomes or analyses appear | Prespecify outcomes and analysis plans |

Safeguards reduce risk; they do not prove bias is absent.

## Confounding Visually

Suppose age affects disease and exposed participants are mostly older.

```text
Exposure --------?--------> Disease
    ^                           ^
    |                           |
    +----------- Age -----------+
```

Age is associated with exposure and independently affects disease. A crude
exposure-disease comparison can therefore mix age with exposure.

## A Numerical Example

Within each age group, exposed and unexposed risks are identical. Because the
exposed group contains mostly older people, the crude comparison looks strong.

```biolang
# Younger: exposed 1/20, unexposed 9/180.
let rr_younger = (1.0 / 20.0) / (9.0 / 180.0)

# Older: exposed 81/180, unexposed 9/20.
let rr_older = (81.0 / 180.0) / (9.0 / 20.0)

# Combined without accounting for age.
let crude_rr = ((1.0 + 81.0) / (20.0 + 180.0)) /
               ((9.0 + 9.0) / (180.0 + 20.0))

{
  younger_risk_ratio: rr_younger,
  older_risk_ratio: rr_older,
  crude_risk_ratio: crude_rr
}
```

Both age-specific risk ratios equal 1, but the crude risk ratio is much larger.
The crude association is explained by the different age composition in this
constructed example.

## How to Address Confounding

Use design before analysis when possible:

- randomization for eligible causal treatment questions;
- restriction or matching when justified;
- careful measurement of plausible confounders;
- stratification, standardization, or regression adjustment;
- sensitivity analysis for unmeasured or imperfectly measured confounding.

Do not adjust blindly for every available variable. A mediator lies on the
causal path; a collider is caused by two variables. Adjusting for either can
answer the wrong question or introduce bias. Draw the assumed causal story
before choosing covariates.

## Association Is Not Causation

Evidence for causality depends on design, time order, alternative explanations,
measurement quality, consistency, magnitude, and subject knowledge. No single
p-value or regression coefficient supplies all of that evidence.

## A Reporting Sentence

> The crude association weakened after age-specific comparison, suggesting
> that age distribution explained much of the unadjusted difference; residual
> confounding remains possible.

> **Do not conclude:** that adjustment proves causality or that an adjusted
> estimate is automatically less biased.

## Quick Check

1. Which problem is usually reduced by increasing an independent sample size?
2. Can perfect precision remove selection bias?
3. Why inspect stratum-specific results before trusting a crude result?

<details>
<summary>Answers</summary>

1. Random sampling uncertainty.
2. No.
3. A third variable may create or hide an association when groups are combined.

</details>

**Next:** [Screening and Diagnostic Tests](simple-09-screening.md). See
[Day 20](day-20.md) for model-based confounder and batch-effect examples.
