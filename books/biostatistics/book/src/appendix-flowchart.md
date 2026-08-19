# Appendix B: Statistical Decision Flowchart

> *The hardest part is translating a scientific question and study design into an analysis. This appendix is a starting map, not an automatic test selector.*

Before choosing a method, write down five things:

1. **Question:** What quantity or contrast do you want to estimate?
2. **Outcome:** Is it continuous, count, binary, ordinal, proportion, or time-to-event?
3. **Experimental unit:** What was independently assigned or sampled—a patient, mouse, culture, library, or cell?
4. **Dependence:** Are observations paired, repeated, clustered, spatial, or time ordered?
5. **Design problems:** Could selection, missingness, confounding, censoring, or measurement error distort the answer?

The tables below narrow the options after those questions are answered. No histogram or normality test can repair a mismatch between the method and the experimental unit.

## The Master Decision Guide

Start here. Find your question type, then follow the table to the right test.

| What are you asking? | Go to section |
|---|---|
| Are two groups different? | [Comparing Two Groups](#comparing-two-groups) |
| Are three or more groups different? | [Comparing Multiple Groups](#comparing-multiple-groups) |
| Are two variables related? | [Associations and Correlations](#associations-and-correlations) |
| Does one variable predict another? | [Regression](#regression) |
| Is there a relationship in categorical data? | [Categorical Data](#categorical-data) |
| How long until an event occurs? | [Time-to-Event Analysis](#time-to-event-analysis) |
| Do I need to reduce dimensionality? | [Dimensionality Reduction](#dimensionality-reduction) |
| Do I need to group similar observations? | [Clustering](#clustering) |

## Comparing Two Groups

Use this when you have one outcome variable and two groups (e.g., control vs. treated, male vs. female, wildtype vs. knockout).

### Step 1: What type is your outcome variable?

| Outcome type | Next step |
|---|---|
| Continuous (expression level, concentration, weight) | Step 2 |
| Counts (number of mutations, colony counts) | Consider a count model; check exposure and overdispersion |
| Binary (alive/dead, present/absent) | See [Categorical Data](#categorical-data) |
| Ordinal (severity scale, ordered category) | Use an ordinal or rank-based method that matches the design |

### Step 2: Are the observations paired or independent?

| Design | Paired? | Example |
|---|---|---|
| Same subjects measured before and after treatment | Yes | Pre/post drug expression |
| Different subjects in each group | No | Treated vs. control mice |
| Matched pairs (e.g., tumor vs. adjacent normal from same patient) | Yes | Tumor/normal tissue pairs |

### Step 3: Choose your test

| Paired? | Working conditions | Variance | Test | BioLang |
|---|---|---|---|---|
| No | Mean difference is the target; no severe influential-point problem | May differ | Welch's t-test | `ttest(a, b, {variance: "welch"})` |
| No | Rank/distribution comparison is appropriate | Not required equal | Mann-Whitney U | `wilcoxon(a, b)` |
| Yes | Mean of paired differences is the target | Applied to differences | Paired t-test | `ttest_paired(a, b)` |
| Yes | Rank-based paired comparison is appropriate | Not required equal | Wilcoxon signed-rank | `wilcoxon_paired(before, after)`; `wilcoxon(a, b)` remains the independent rank-sum test |

> **Key insight:** Use `{variance: "welch"}` for Welch's independent-samples t-test; the two-argument form remains pooled for backward compatibility. A rank test is not simply a “t-test without normality”; it changes what feature of the distributions is being compared and needs its own interpretation.

### How to check normality

```bio
import "statistics" as stat

let data = [2.3, 4.1, 3.7, 5.2, 4.8, 3.1, 6.0, 4.4]
let report = stat.explore(data)
stat.distribution_plot(data, {title: "Distribution check"})
```

For an independent-group mean comparison, inspect each group and the fitted model. For a paired t-test, inspect the **paired differences**, because those are what the test analyses. Use `normal_qq_plot(values)` as one visual diagnostic. The older `qq_plot(p_values)` remains specifically for genomic p-value Q-Q plots.

> **Common pitfall:** A normality-test p-value is not a traffic light. With small samples it may miss important departures; with very large samples it may detect harmless ones. Combine Q-Q plots, raw points, influence checks, study design, and sensitivity analysis.

## Comparing Multiple Groups

Use this when you have three or more groups (e.g., three drug doses, four tissue types, five time points).

| Normal? | Equal variance? | Design | Test | BioLang |
|---|---|---|---|---|
| Yes | Yes | Independent groups | One-way ANOVA | `anova(groups)` |
| Yes | No | Independent groups | Welch's ANOVA | `anova(groups, {variance: "welch"})` |
| No | — | Independent groups | Kruskal-Wallis rank/distribution comparison | `kruskal_wallis(groups)` |
| Yes | — | Repeated measures | Repeated-measures ANOVA | Not currently a builtin |
| No | — | Repeated measures | Friedman test | Not currently a builtin |
| Yes | — | Two factors | Two-way ANOVA | Use an explicitly specified regression model; `anova(groups)` is one-way only |

### Post-hoc Tests

When ANOVA is significant, you know *some* groups differ but not *which* ones. Use post-hoc tests:

| Test | When to use | BioLang |
|---|---|---|
| Tukey HSD | All pairs after a defensible classical equal-variance ANOVA | `tukey_hsd(groups)` |
| Dunnett | Compare all groups to a single control | Not currently a builtin; adjusted pairwise tests are not Dunnett's joint procedure |
| Dunn test | Post-hoc for Kruskal-Wallis | Not currently a builtin; separate rank-sum tests are not Dunn's pooled-rank procedure |
| Adjusted pairwise mean tests | Pair-specific Welch or pooled comparisons | `pairwise_ttest(groups, {variance: "welch", adjust: "holm"})` |

> **Key insight:** ANOVA is an *omnibus* test—it addresses whether the group means are all compatible with equality under the model. If the scientific question concerns particular groups, report pre-planned contrasts or multiplicity-aware follow-up comparisons with effect estimates and intervals.

## Associations and Correlations

Use this when you have two continuous variables and want to know if they are related (e.g., gene expression vs. methylation, age vs. telomere length).

| Data characteristics | Test | BioLang |
|---|---|---|
| Both variables roughly normal, linear relationship | Pearson correlation | `cor(x, y)` |
| Non-normal or ordinal data, monotonic relationship | Spearman correlation | `spearman(x, y)` |
| Ordinal data with ties | Kendall tau | `kendall(x, y)` |
| Partial correlation (controlling for a third variable) | Partial correlation | `cor(x, y)` after residualizing on z |

### Interpreting Correlation Magnitude

The following labels are rough descriptions, not biological decision thresholds. A correlation of 0.2 may matter for a widespread exposure, while 0.8 may still be useless for individual prediction. Always show the scatter plot, sample size, interval, and scientific units.

| |r| value | Interpretation |
|---|---|
| 0.0 - 0.1 | Negligible |
| 0.1 - 0.3 | Weak |
| 0.3 - 0.5 | Moderate |
| 0.5 - 0.7 | Strong |
| 0.7 - 1.0 | Very strong |

> **Common pitfall:** Correlation does not imply causation, but more subtly, *absence* of Pearson correlation does not imply absence of relationship. Pearson only detects linear associations. Two variables can have a perfect quadratic relationship with r = 0. Always plot your data.

## Categorical Data

Use this when both your variables are categorical (e.g., mutation status vs. disease outcome, genotype vs. phenotype).

| Design | Expected cell counts | Test | BioLang |
|---|---|---|---|
| 2x2 table, large samples | All expected >= 5 | Chi-square test of independence | Constructing expected contingency counts is not currently a dedicated builtin |
| 2x2 table, small samples | Any expected < 5 | Fisher's exact test | `fisher_exact(a, b, c, d)` |
| Larger than 2x2 | All expected >= 5 | Chi-square test of independence | Not currently a dedicated builtin; `chi_square()` accepts flat observed/expected lists for goodness-of-fit |
| Larger than 2x2, small samples | Any expected < 5 | Fisher-Freeman-Halton | Not currently a builtin; `fisher_exact()` is 2x2 only |
| Paired categorical data | — | McNemar's test | Not currently a builtin |
| Trend across ordered categories | — | Cochran-Armitage trend test | Not currently a builtin |

### Measures of Association for Categorical Data

| Measure | Use case | BioLang |
|---|---|---|
| Odds ratio | 2x2 tables, case-control studies | `(a*d) / (b*c)` (inline) |
| Relative risk | 2x2 tables, cohort studies | `(a/(a+b)) / (c/(c+d))` (inline) |
| Cramer's V | Any size contingency table | Compute from chi-square statistic |

## Regression

Use this when you want to predict an outcome from one or more predictor variables.

| Outcome type | Number of predictors | Test | BioLang |
|---|---|---|---|
| Continuous | 1 | Simple linear regression | `lm(x, y)` |
| Continuous | Multiple | Multiple linear regression | `lm(~y ~ x1 + x2 + x3, table)` |
| Binary (0/1) | Any | Logistic regression | `glm(~y ~ x, table, "binomial")` |
| Count | Any | Poisson regression | `glm(~y ~ x, table, "poisson")` |
| Count, overdispersed | Any | External negative-binomial tool or package | Not a current `glm` family |
| Continuous, clustered data | Any | Stratified models or an external mixed-effects tool | Fit and report per group |

### Checking Regression Assumptions

```bio
let age = [34, 41, 52, 59, 63, 70]
let expression = [8.1, 8.8, 9.4, 10.2, 10.6, 11.1]
let model = lm(age, expression)

# Compute residuals from the fitted line
let residuals = range(0, len(age))
  |> map(|i| expression[i] - (model.intercept + model.slope * age[i]))
histogram(residuals, {title: "Residual Distribution"})
print("R-squared: " + str(round(model.r_squared, 3)))
```

> **Common pitfall:** Adding more predictors always improves R-squared, even if the predictors are noise. Use adjusted R-squared or AIC/BIC for model comparison. Report both R-squared and adjusted R-squared.

## Time-to-Event Analysis

Use this when your outcome is the time until something happens (death, relapse, response) and some observations are censored (the event has not yet occurred).

| Question | Method | BioLang |
|---|---|---|
| Estimate survival curve | Kaplan-Meier | Sort event times, compute stepwise survival |
| Compare survival between two groups | Log-rank test | `log_rank_test(times_a, events_a, times_b, events_b)` |
| Compare survival, multiple groups | Multi-group log-rank test | Not currently a builtin |
| Adjust for covariates | Cox proportional hazards | `cox_ph(time, event, covariates)` |
| Estimate median survival | From sorted times | `sort(times)[len(times) / 2]` |

> **Clinical relevance:** Some clinical trials use a Cox-model hazard ratio as an efficacy measure. A hazard ratio of 0.65 describes a modeled instantaneous event-rate ratio, conditional on the model; it is not automatically a 35% reduction in an individual's probability of experiencing the event. Report the interval, absolute survival summaries, censoring information, and proportional-hazards assessment.

## Dimensionality Reduction

Use this when you have many variables (genes, proteins, metabolites) and want to find the main patterns.

| Goal | Method | BioLang |
|---|---|---|
| Find linear combinations that maximize variance | PCA | `pca(data)` |
| Visualize PCA results | PCA plot | `pca_plot(result, {title: "PCA"})` |

> **Key insight:** A fixed PCA implementation and input generally reproduces the same subspace, although component signs and nearly tied components can differ. t-SNE and UMAP commonly use stochastic or approximate steps, so record seeds, parameters, software versions, and the input representation.

## Clustering

Use this when you want to group similar observations (samples, genes, cells) together.

| What you know | Method | BioLang |
|---|---|---|
| Number of clusters (k) | k-means | `kmeans(data, 3)` |
| Want a hierarchy of clusters | Hierarchical clustering | `hclust(data, "ward")` |
| Irregular cluster shapes | DBSCAN | `dbscan(data, 0.5, 5)` |
| Want to estimate k | Silhouette / Elbow | Loop over k, compute `kmeans(data, k).silhouette` |

## Multiple Testing Correction

Use this whenever you perform more than one statistical test on the same dataset.

| Method | Controls | Strictness | BioLang |
|---|---|---|---|
| Bonferroni | Family-wise error rate | Most conservative | `p_adjust(pvals, "bonferroni")` |
| Holm | Family-wise error rate | Less conservative | `p_adjust(pvals, "holm")` |
| Benjamini-Hochberg | False discovery rate | Moderate | `p_adjust(pvals, "BH")` |
| Benjamini-Yekutieli | FDR under dependence | Conservative FDR | Not currently a builtin |
| Permutation | Empirical null | Gold standard | Inline loop with `shuffle()` |

> **Key insight:** In exploratory genomics, Benjamini-Hochberg FDR control is common because it balances discovery and false-discovery control. Family-wise methods such as Holm or Bonferroni answer a stricter question and may be appropriate when even one false positive is costly. Choose the error criterion from the study goal rather than habit.

## Quick Reference: Common Biological Scenarios

| Scenario | Recommended test | BioLang |
|---|---|---|
| Gene expression, treated vs. control | Welch's t-test | `ttest(treated, control, {variance: "welch"})` |
| Gene expression across 4 tissues | One-way ANOVA | `anova([tissue1, tissue2, tissue3, tissue4])` |
| Mutation frequency in cases vs. controls | Fisher's exact test | `fisher_exact(a, b, c, d)` |
| Survival by treatment arm | Log-rank comparison with censoring | `log_rank_test(times_a, events_a, times_b, events_b)` |
| 20,000 gene differential expression | t-test + BH correction | `p_adjust(pvals, "BH")` |
| Sample clustering from RNA-seq | PCA + hierarchical clustering | `pca(data)` then `hclust(scores)` |
| Correlation: expression vs. methylation | Spearman (often non-linear) | `spearman(expr, meth)` |
| GWAS: genotype vs. phenotype | Logistic regression + BH | `glm(~pheno ~ geno, tbl, "binomial")` |
| Clinical outcome predictors | Regression model | `lm(~outcome ~ age + stage + treatment, tbl)` |
| Sample size for planned experiment | Power analysis | Compute with `qnorm()` and effect size |
