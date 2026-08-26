# Start Simple 1: Centre - Mean, Median, and Mode

> **In one sentence:** Centre is a chosen way to describe what is typical; the
> mean, median, and mode answer different versions of that question.

## The Question

Five samples have protein concentrations of 4, 5, 5, 6, and 30 ng/mL. What
single value should describe a typical sample?

```text
4    5    5    6                              30
|----|----|----|-------------------------------|
          ^              ^
       median 5        mean 10
       mode 5
```

The isolated value 30 pulls the mean to the right. It barely changes the
median, and it does not change the most frequent value.

## Three Different Centres

| Centre | Plain meaning | Best question | Sensitive to extremes? |
|---|---|---|---|
| **Mean** | Add every value and divide equally | What is the arithmetic balance point or total per observation? | Yes |
| **Median** | Sort values and take the halfway value | What does a typical ranked observation look like? | Much less |
| **Mode** | The most frequent category or value | What occurs most often? | Usually no |

There is no universally best centre. Choose it from the scientific question,
then inspect the values or a plot to see what that single number hides.

## Calculate It in BioLang

```biolang
let protein_values = [4.0, 5.0, 5.0, 6.0, 30.0]

let centre = {
  mean: mean(protein_values),
  median: median(protein_values),
  overview: summary(protein_values)
}

centre
```

`protein_values` is deliberately named so it cannot be confused with a
function such as `protein()`. `summary()` keeps the minimum and maximum beside
the centres, making the distant value harder to overlook.

For this tiny dataset, the mode is visibly 5 because it occurs twice. Do not
force a mode from continuous measurements when nearly every value is unique.

## Which One Should You Report?

- Use **mean with SD** when the arithmetic balance point is meaningful and a
  few values are not dominating the story.
- Use **median with IQR** when the halfway observation is more meaningful, the
  data are skewed, or valid extreme values make the mean misleading.
- Use **mode with counts or percentages** for categories such as blood group.
- Always keep a dot plot, histogram, or the raw values nearby.

### Other Kinds of Mean

These are not interchangeable decorations:

| Mean | Use it when |
|---|---|
| Arithmetic mean | Values combine by addition |
| Geometric mean | Positive values combine multiplicatively, such as growth factors or fold changes |
| Weighted mean | Observations have justified, known contributions or sampling weights |
| Harmonic mean | Averaging rates over equal amounts of work, such as equal-distance speeds |

If you cannot explain what the weights or multiplicative scale mean, use the
raw data and seek advice before choosing a specialized mean.

## A Reporting Sentence

> The median protein concentration was 5 ng/mL; the arithmetic mean was 10
> ng/mL because one sample measured 30 ng/mL.

That sentence reports the value and explains why two centres disagree.

> **Do not conclude:** that 30 is an error merely because it is far away. Check
> the sample, instrument, units, and protocol before excluding any observation.

## Quick Check

1. Which centre divides ordered observations into two equal halves?
2. Which centre preserves the arithmetic total, because `mean x count = sum`?
3. What does a large mean-median difference ask you to inspect?

<details>
<summary>Answers</summary>

1. The median.
2. The arithmetic mean.
3. Skewness, extreme values, mixed groups, data errors, or another unusual
   shape. The difference is a clue, not a diagnosis.

</details>

**Next:** [Spread, Shape, Outliers, and Logs](simple-02-spread-shape.md). For
more depth, see [Day 2](day-02.md) and [Day 3](day-03.md).
