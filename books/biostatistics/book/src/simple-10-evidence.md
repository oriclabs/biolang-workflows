# Start Simple 10: Read the Evidence

> **In one sentence:** A trustworthy conclusion depends on the question,
> design, measurement, effect size, uncertainty, bias, and applicability - not
> on one impressive number.

## Turn the Claim into a Question

For an intervention question, PICO is a useful starting frame:

| Part | Ask |
|---|---|
| Population | Who is the decision about? |
| Intervention or exposure | What is being done or observed? |
| Comparator | Compared with what? |
| Outcome and time | What matters, and when is it measured? |

"Does treatment X work?" is too vague. "Among adults meeting these eligibility
criteria, does X compared with usual care reduce hospitalization within 90
days?" can guide design, data, and interpretation.

## Read a Result in Layers

```text
1. What was the question?
2. Who was studied and who was not?
3. How were groups formed?
4. What was measured, and was measurement fair?
5. How large was the effect in absolute and relative terms?
6. How uncertain was it?
7. What bias, confounding, missingness, or selective reporting could remain?
8. Does it apply to the decision in front of me?
```

## Estimate, Interval, and Null Value

Suppose a study reports a risk ratio of 0.78 with a 95% confidence interval from
0.64 to 0.95. A ratio null value is 1.

```biolang
let result = {
  estimate: 0.78,
  confidence_low: 0.64,
  confidence_high: 0.95,
  null_value: 1.0
}

let interval_crosses_null =
  result.confidence_low <= result.null_value and
  result.confidence_high >= result.null_value

{
  relative_reduction: 1.0 - result.estimate,
  interval_crosses_null: interval_crosses_null
}
```

The point estimate suggests a 22% relative reduction, and this interval does
not include 1. But you still need baseline risks to understand the absolute
benefit. Reducing risk from 50% to 39% and from 0.5% to 0.39% gives the same
relative ratio but very different practical consequences.

## Internal and External Validity

**Internal validity** asks whether the comparison supports the claimed effect
for the people studied. Consider randomization or confounding, allocation,
masking, measurement, missing outcomes, adherence, analysis population, and
selective reporting.

**External validity** asks whether the evidence applies elsewhere. Compare the
study population, setting, intervention, comparator, outcome, follow-up, and
available care with the real decision.

A study can be internally strong but poorly applicable to another population.
It can also be highly representative but too biased internally to support its
claim.

## Evidence Synthesis Helps, but Cannot Repair Weak Inputs

A systematic review uses an explicit search and appraisal process. A
meta-analysis may combine sufficiently compatible estimates. A forest plot
shows each estimate and interval plus any pooled result.

Pooling increases precision only when the studies, outcomes, and models can be
meaningfully combined. It does not erase publication bias, confounding, poor
measurement, duplicate populations, or incompatible questions.

## A Minimal Appraisal Note

Write one sentence for each:

```text
Question:       the exact population, comparison, outcome, and time
Design:         how groups and measurements were created
Result:         absolute effect, relative effect, and interval
Limitations:    the most credible threats to the conclusion
Applicability:  how closely the study matches the intended decision
Decision:       what the evidence supports now, and what remains uncertain
```

## A Reporting Sentence

> The estimated risk ratio favored treatment, but its practical importance
> depends on baseline absolute risk; interpretation also considers study
> conduct, missing outcomes, harms, and applicability to the target population.

> **Do not conclude:** that a randomized trial is automatically flawless, that
> an observational study is automatically useless, or that a meta-analysis is
> automatically the highest-quality answer.

## Quick Check

1. Which validity asks whether the comparison is credible within the study?
2. Why report absolute risk alongside a relative effect?
3. Can a precise pooled estimate correct biased source studies?

<details>
<summary>Answers</summary>

1. Internal validity.
2. The same relative effect can imply very different numbers helped or harmed
   at different baseline risks.
3. No. Greater numerical precision does not remove systematic error.

</details>

## Where to Continue

You have completed Start Simple. Choose the full chapters that match your work:

| Need | Continue with |
|---|---|
| Describe and visualize data | [Days 2-3](day-02.md) |
| Confidence intervals and tests | [Days 6-12](day-06.md) |
| Regression and binary outcomes | [Days 13-16](day-13.md) |
| Survival and study planning | [Days 17-19](day-17.md) |
| Confounding and batch effects | [Day 20](day-20.md) |
| Omics and high-dimensional data | [Days 21-30](day-21.md) |
| Choose a method quickly | [Decision Flowchart](appendix-flowchart.md) |

Or return to [Choose Your Path](start-simple.md).
