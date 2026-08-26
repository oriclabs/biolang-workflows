# Start Simple 9: Screening and Diagnostic Tests

> **In one sentence:** Sensitivity and specificity describe test behavior;
> predictive values describe what a result means in a particular population.

## Screening Is Not Diagnosis

Screening looks for possible disease in people without recognized symptoms.
It is the beginning of a pathway, not a final diagnosis. A useful screening
programme needs more than an accurate test: an important condition, a useful
detectable stage, acceptable follow-up, effective action, fair access, and a
favorable balance of benefit, harm, and cost.

## Start with the 2x2 Test Table

| | Disease present | Disease absent |
|---|---:|---:|
| Test positive | True positive | False positive |
| Test negative | False negative | True negative |

```text
Sensitivity = true positives / everyone with disease
Specificity = true negatives / everyone without disease
PPV         = true positives / everyone who tests positive
NPV         = true negatives / everyone who tests negative
```

Sensitivity and specificity condition on disease status. PPV and NPV condition
on the observed test result.

## Why Prevalence Changes Meaning

Consider 10,000 screened people, disease prevalence 1%, sensitivity 90%, and
specificity 95%.

```biolang
let population = 10000.0
let prevalence = 0.01
let sensitivity = 0.90
let specificity = 0.95

let disease_present = population * prevalence
let disease_absent = population - disease_present
let true_positive = disease_present * sensitivity
let false_negative = disease_present - true_positive
let true_negative = disease_absent * specificity
let false_positive = disease_absent - true_negative

let ppv = true_positive / (true_positive + false_positive)
let npv = true_negative / (true_negative + false_negative)

{
  true_positive: true_positive,
  false_positive: false_positive,
  false_negative: false_negative,
  true_negative: true_negative,
  positive_predictive_value: ppv,
  negative_predictive_value: npv
}
```

Even a fairly specific test can produce more false positives than true
positives when the condition is uncommon. That is arithmetic, not test failure.
Testing a different-risk population changes PPV and NPV even if sensitivity
and specificity stay the same.

## Thresholds Trade Errors

For many continuous tests, changing the positive threshold changes both
sensitivity and specificity.

```text
lower threshold  -> more positives -> higher sensitivity, lower specificity
higher threshold -> fewer positives -> lower sensitivity, higher specificity
```

The sensible threshold depends on consequences. Missing a treatable dangerous
condition and causing an unnecessary follow-up procedure have different harms.
A ROC curve summarizes threshold trade-offs but does not choose the threshold
or establish clinical utility.

## Spectrum and Verification Matter

Test performance can change with disease stage, age, comorbidity, sample
quality, setting, and how the reference diagnosis is established. If only
positive tests receive definitive verification, accuracy estimates can be
biased.

## A Reporting Sentence

> In this 1%-prevalence screening scenario, predictive values were calculated
> from sensitivity, specificity, and the target population; a positive screen
> therefore required confirmatory assessment rather than being treated as a
> diagnosis.

> **Do not conclude:** that a test with 90% sensitivity means a person with a
> positive result has a 90% probability of disease.

## Quick Check

1. Which quantity answers: among positive tests, how many truly have disease?
2. What usually happens to PPV when prevalence decreases?
3. Does a high area under the ROC curve prove that screening improves health?

<details>
<summary>Answers</summary>

1. Positive predictive value.
2. PPV usually decreases, with sensitivity and specificity held constant.
3. No. Benefits, harms, follow-up, treatment, implementation, and the target
   population must also be evaluated.

</details>

**Next:** [Reading Health Evidence](simple-10-evidence.md). Logistic modelling
and ROC curves are covered in [Day 16](day-16.md).
