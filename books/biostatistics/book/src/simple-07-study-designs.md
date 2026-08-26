# Start Simple 7: Sampling and Study Designs

> **In one sentence:** The study design determines which comparisons are
> possible, when measurements occur, and how far a result can be generalized.

## Sampling and Design Are Different Decisions

**Sampling** asks who enters the study. **Study design** asks what is measured,
when it is measured, and whether an exposure or treatment is assigned.

```text
Target population
       |
       | sampling
       v
Study participants
       |
       | study design
       v
Measurements and comparisons
```

A randomized treatment assignment does not guarantee a representative sample.
A representative survey does not create a randomized experiment.

## Common Sampling Approaches

| Approach | How units enter | Main caution |
|---|---|---|
| Simple random | Draw from a complete sampling frame | Requires a usable frame |
| Stratified | Draw within predefined groups | Analysis should respect the design and weights |
| Systematic | Random start, then every kth unit | List periodicity can create bias |
| Cluster | Draw groups such as clinics or schools | Units within a cluster are correlated |
| Convenience | Recruit what is readily available | Unknown selection mechanism limits generalization |

BioLang can demonstrate a random draw, but software cannot make an incomplete
sampling frame representative:

```biolang
let sampling_frame = [101, 102, 103, 104, 105, 106, 107, 108,
                      109, 110, 111, 112, 113, 114, 115, 116]
let selected_ids = sample(sampling_frame, 6)
selected_ids
```

Record the seed for reproducible analysis. For real recruitment or treatment
allocation, use an approved procedure with allocation safeguards rather than
an informal script.

## Four Basic Study Designs

| Design | Starting point | Time direction | Useful for | Important limitation |
|---|---|---|---|---|
| Cross-sectional | A sample now | One time point | Prevalence and current associations | Temporality is often unclear |
| Cohort | Exposure or eligibility | Forward through follow-up, sometimes using existing records | Incidence, risks, prognosis | Attrition and confounding |
| Case-control | Outcome status | Look back for prior exposure | Rare outcomes or long latency | Recall and selection bias; sampled risks unavailable |
| Randomized trial | Eligible participants | Assign treatment, then follow | Estimating treatment effects under the trial conditions | Non-adherence, missing outcomes, ethics, and generalizability |

```text
Cross-sectional:   exposure? ---- outcome?     measured now
Cohort:            exposure  ----------------> outcome
Case-control:      exposure? <---------------- case/control
Randomized trial:  randomize ----------------> outcome
```

## Match the Claim to the Design

- A cross-sectional association usually cannot show which came first.
- A cohort establishes time order more clearly but can remain confounded.
- A case-control study efficiently studies rare outcomes but usually estimates
  an odds ratio rather than population risk from sampled counts.
- Randomization helps balance causes on average, but execution, adherence,
  missingness, and analysis still matter.

## Questionnaires Are Measurements

A questionnaire is not merely a list of questions. Define the construct,
choose a reference period, avoid leading or double questions, pilot the wording,
and assess reliability and validity. Missing responses and non-response can be
systematic rather than harmless.

```text
Bad:  "How often do poor diet and inactivity make you tired?"
      two exposures + judgment + one outcome

Better: separate diet, activity, and fatigue into defined questions
        with a stated reference period.
```

## A Reporting Sentence

> This was a cohort study of a defined sampled population; exposure preceded
> outcome measurement, but the observational comparison remains vulnerable to
> confounding and loss to follow-up.

> **Do not conclude:** that "prospective" means randomized, or that
> randomization makes every later analysis unbiased.

## Quick Check

1. Which design begins by selecting cases and controls?
2. Which design is suited to estimating point prevalence?
3. Does a random sample automatically imply randomized treatment assignment?

<details>
<summary>Answers</summary>

1. Case-control.
2. Cross-sectional.
3. No. Sampling and treatment assignment are separate processes.

</details>

**Next:** [Bias, Confounding, and Fair Comparisons](simple-08-bias-confounding.md).
For fuller design and power discussion, see [Day 18](day-18.md).
