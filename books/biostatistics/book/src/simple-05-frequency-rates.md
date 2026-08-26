# Start Simple 5: Counts, Proportions, and Disease Rates

> **In one sentence:** A health number becomes interpretable only when its
> denominator, population, place, and time period are clear.

## Four Related Numbers

Imagine a town of 10,000 people. At the start of a year, 200 people have a
condition. During the year, 40 new cases occur among people who could develop
it.

| Measure | Plain question | Basic form |
|---|---|---|
| Count | How many? | number of events or people |
| Proportion | What fraction of the whole? | part / whole |
| Prevalence | Who has the condition at a stated time or during a stated period? | existing cases / population |
| Incidence | Who developed the condition while at risk? | new cases / population at risk, or person-time at risk |

```text
Prevalence is a photograph:       who has it now?
Incidence is a video:             who developed it over time?

Existing at start:  O O O O
New during follow-up:       + +
```

## Calculate the Measures

```biolang
let population = 10000.0
let existing_cases = 200.0
let new_cases = 40.0
let person_years_at_risk = 9500.0

let prevalence = existing_cases / population
let incidence_rate = new_cases / person_years_at_risk

{
  prevalence_percent: prevalence * 100.0,
  incidence_per_1000_person_years: incidence_rate * 1000.0
}
```

The prevalence is expressed as a percentage of a defined population at a
defined time. The incidence rate is expressed per 1,000 person-years, not per
1,000 people with no time attached.

## Risk and Rate Are Not Identical

**Risk** is a probability over a stated interval. It has a denominator of
people initially at risk.

**Rate** allows different follow-up times. Its denominator is person-time, such
as person-years at risk. A rate can exceed 1 per person-year when recurrent
events are counted; a probability cannot exceed 1.

Always ask:

```text
Numerator:    what event was counted?
Denominator:  who or what could contribute?
Time:         at what point or over what interval?
Population:   which people, samples, or organisms?
```

## Crude and Specific Rates

A crude rate combines everyone. A specific rate describes a subgroup, such as
an age group. Two regions can have different crude death rates simply because
one has an older population.

Age standardization answers a comparison question: what rates would we expect
if both populations had the same age distribution?

```biolang
# Two age-specific rates, combined using a common reference population.
let rate_younger = 0.002
let rate_older = 0.030
let reference_weight_younger = 0.70
let reference_weight_older = 0.30

let standardized_rate =
  rate_younger * reference_weight_younger +
  rate_older * reference_weight_older

standardized_rate * 1000.0
```

This returns a standardized rate per 1,000 under the chosen reference weights.
It is a comparison tool, not the observed crude rate of either population.

## Why Prevalence Can Be High Without High Incidence

Prevalence is influenced by both new cases and duration:

```text
many new cases + long duration   -> often high prevalence
few new cases  + long duration   -> prevalence may still be high
many new cases + rapid recovery  -> prevalence may remain modest
```

Therefore prevalence alone does not tell you how quickly new disease is
appearing.

## A Reporting Sentence

> At the stated date, prevalence was 2.0% in the defined population; during the
> year, the incidence rate was reported per 1,000 person-years at risk.

> **Do not conclude:** that two rates are comparable until their case
> definitions, denominators, observation times, and population structures are
> comparable.

## Quick Check

1. Which measure counts new cases over follow-up?
2. Why is "40 cases" incomplete as a rate?
3. Why might age standardization change a regional comparison?

<details>
<summary>Answers</summary>

1. Incidence risk or incidence rate, depending on the denominator.
2. It lacks the population or person-time that could produce those cases.
3. It removes differences caused merely by the regions having different age
   compositions under the chosen reference population.

</details>

**Next:** [One 2x2 Table, Many Useful Measures](simple-06-two-by-two.md).
