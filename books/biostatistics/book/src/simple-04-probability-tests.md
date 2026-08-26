# Start Simple 4: Probability, Tests, and p-Values

> **In one sentence:** A p-value measures compatibility with a stated null
> model; it is not the probability that the null hypothesis is true.

## Begin with the Question

Suppose treated and control samples have different mean expression. Three
questions must stay separate:

1. **Estimate:** how large is the observed difference?
2. **Uncertainty:** how precisely was that difference estimated?
3. **Test:** would a difference at least this extreme be surprising under a
   specified no-difference model?

The p-value answers only the third question.

```text
Observed difference
        |
        v
Assume the null model for comparison
        |
        v
How often would results this extreme or more extreme occur?
        |
        v
      p-value
```

## Probability Is Conditional

`P(data this extreme | null model)` is not the same as
`P(null model | observed data)`. Reversing the condition changes the question.

A small p-value says the data are difficult to reconcile with the tested null
model and its assumptions. It does not by itself say the effect is large,
important, causal, correctly measured, or likely to reproduce.

## A Small BioLang Comparison

```biolang
let control_values = [7.8, 8.1, 8.4, 7.9, 8.3, 8.0, 8.2, 7.7]
let treated_values = [8.6, 8.9, 8.5, 9.1, 8.8, 8.7, 9.0, 8.4]

let difference = mean(treated_values) - mean(control_values)
let test = ttest(treated_values, control_values, {variance: "welch"})

{
  mean_difference: difference,
  p_value: test.p_value,
  test_statistic: test.statistic
}
```

Welch's test compares independent group means without requiring equal group
variances. It still assumes the experimental units are independent and that a
difference in means is the intended question.

## What Does 0.05 Mean?

`0.05` is a convention, not a law of nature.

```text
p = 0.049     and     p = 0.051
       nearly the same evidence
       not opposite scientific truths
```

Predeclare an error threshold when a formal decision is required, but report
the exact p-value with the effect estimate and confidence interval. Scientific
importance comes from the size, uncertainty, consequences, design, and prior
evidence - not from which side of 0.05 a result lands.

## False Alarms and Missed Effects

| Reality and decision | Meaning |
|---|---|
| No relevant effect, but claim one | Type I error or false alarm |
| Relevant effect, but fail to detect it | Type II error or missed effect |
| Detect a real effect | Power is the long-run probability of doing this under a specified alternative |

Failure to reject a null hypothesis is not proof of no effect. A wide interval
may still include both meaningful benefit and meaningful harm.

## Multiple Testing

If 10,000 genes are tested at an unadjusted 0.05 threshold under complete null
conditions, many small p-values are expected by chance. Use a planned
multiple-testing method such as false discovery rate control and validate
important findings independently.

## A Reporting Sentence

> The treated-control mean difference was reported with its confidence
> interval; Welch's test supplied a p-value describing compatibility with the
> no-mean-difference model.

> **Do not conclude:** "there is no difference" merely because `p > 0.05`, or
> "the treatment works" merely because `p < 0.05`.

## Quick Check

1. Is a p-value the probability that the null hypothesis is true?
2. Can a tiny, unimportant effect have a small p-value in a very large study?
3. What three items should accompany a test result?

<details>
<summary>Answers</summary>

1. No.
2. Yes.
3. The effect estimate, its uncertainty or confidence interval, and the study
   design and assumptions. Report the p-value as supporting information.

</details>

**Next:** [Counts, Proportions, and Disease Rates](simple-05-frequency-rates.md).
For more depth, see [Day 7](day-07.md), [Day 8](day-08.md), and
[Day 12](day-12.md).
