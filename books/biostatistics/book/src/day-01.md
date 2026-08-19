# Day 1: Why Statistics? The Story Your Data Is Trying to Tell

> **Start here**
> - **In one sentence:** Statistics helps us learn from incomplete, variable data without pretending that uncertainty has disappeared.
> - **Look for:** the biological question, the experimental unit, possible bias, the size of the effect, and its uncertainty.
> - **Use this when:** reading a paper, planning an experiment, checking data, or making a claim from observations.
> - **Do not conclude:** that a calculation can rescue a poorly designed study or turn association into causation.

<div class="day-meta">
<span class="badge">Day 1 of 30</span>
<span class="badge">No Prerequisites</span>
<span class="badge">~45 min reading</span>
<span class="badge">Motivation & Context</span>
</div>

## The Problem

Suppose a small early study reports that a treatment helped. The estimated benefit looks large, but the confidence interval is wide. A later study obtains a smaller estimate.

Neither result has to be fraudulent or incompetent. Small samples naturally produce estimates that move around more from study to study. Bias, different patient populations, measurement choices, and analysis decisions can add further disagreement.

Statistics helps you ask four practical questions: **What was measured? How variable is it? How uncertain is the estimate? Could the study design have distorted it?** A p-value is only one part of those answers.

Over 30 days, this book develops that reasoning with biological examples and BioLang code. The aim is not to memorize test names. It is to understand what a result can—and cannot—support.

## What Is Statistics?

Statistics is a way to learn from data while keeping uncertainty visible. It helps us summarize observations, compare explanations, estimate effects, and describe how much an answer may change with another sample.

Biological data vary for many reasons. Organisms differ, instruments add measurement error, sampling captures only part of a population, and experiments may contain bias or confounding. Statistics cannot erase those problems. It can help us see them and avoid claiming more than the study supports.

## Why Results May Not Reproduce

Two studies can disagree because of ordinary sampling variation, different protocols or populations, weak measurement, selective reporting, confounding, or an analysis that does not match the design. Replication is therefore a new measurement of the claim—not a ceremonial repeat of the same calculation.

Statistical reasoning helps by making the question, experimental unit, analysis choices, uncertainty, and limitations explicit. It cannot compensate for missing controls or biased sampling, but it makes those weaknesses harder to hide.

### Trying Many Analyses and Reporting One

"P-hacking" means allowing the observed result to guide repeated analysis choices, then reporting only the analysis that crossed a desired threshold. A researcher might:

- Test 15 different subgroups and report only the one with p < 0.05
- Remove outliers selectively until the result becomes significant
- Try multiple statistical tests and report whichever gives the smallest p-value
- Add or remove covariates until the "right" answer appears
- Decide when to stop collecting data based on whether the current result is significant

Each choice may sound defensible in isolation. Together, undisclosed flexibility increases the chance of finding an apparently persuasive result even when no repeatable effect exists. Pre-specification, transparent sensitivity analyses, and reporting all relevant outcomes make the evidence easier to evaluate.

### Underpowered Studies

Small studies usually estimate effects imprecisely. Their power depends on the effect size, biological variation, design, measurement error, and chosen decision rule—not on sample size alone. When only threshold-crossing results are published, the reported effects from small studies can also look larger than the underlying effects.

### The Garden of Forking Paths

Researchers make dozens of analytical decisions: how to clean the data, which variables to include, how to handle missing values, which test to use, whether to transform the data, how to define the outcome. Each decision is a fork in the path, and different choices lead to different results. When these choices are made after seeing the data (rather than pre-specified in an analysis plan), the researcher unconsciously navigates toward significance.

These problems are usually workflow problems, not evidence of deliberate misconduct. Clear study questions, recorded decisions, appropriate models, readable diagnostics, and reproducible code all help.

> **Key insight:** A result is easier to trust when another person can see what was planned, what changed, what was measured, and how sensitive the conclusion is to reasonable choices.

## Signal vs. Noise

Here is the most fundamental question in statistics: **Is the pattern I see real, or could it have happened by chance?**

<div style="text-align: center; margin: 2em 0;">
<svg width="680" height="320" viewBox="0 0 680 320" xmlns="http://www.w3.org/2000/svg" style="background: #fafbfc; border: 1px solid #e5e7eb; border-radius: 8px;">
  <defs>
    <marker id="arrow1" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><path d="M0,0 L8,3 L0,6 Z" fill="#6b7280"/></marker>
  </defs>
  <text x="340" y="28" text-anchor="middle" font-size="15" font-weight="bold" fill="#1e293b">Signal vs. Noise</text>
  <!-- Left panel: noisy data with no real trend -->
  <rect x="30" y="42" width="290" height="240" rx="6" fill="#f8fafc" stroke="#e2e8f0"/>
  <text x="175" y="62" text-anchor="middle" font-size="12" font-weight="600" fill="#dc2626">Small Sample (n=10) — Noise Mimics Signal</text>
  <!-- Axes -->
  <line x1="60" y1="250" x2="295" y2="250" stroke="#9ca3af" stroke-width="1.5"/>
  <line x1="60" y1="250" x2="60" y2="75" stroke="#9ca3af" stroke-width="1.5"/>
  <text x="175" y="272" text-anchor="middle" font-size="11" fill="#6b7280">Measurement</text>
  <text x="45" y="165" text-anchor="middle" font-size="11" fill="#6b7280" transform="rotate(-90,45,165)">Response</text>
  <!-- Scattered noisy points -->
  <circle cx="80" cy="200" r="5" fill="#ef4444" opacity="0.7"/><circle cx="105" cy="120" r="5" fill="#ef4444" opacity="0.7"/>
  <circle cx="125" cy="230" r="5" fill="#ef4444" opacity="0.7"/><circle cx="150" cy="105" r="5" fill="#ef4444" opacity="0.7"/>
  <circle cx="170" cy="180" r="5" fill="#ef4444" opacity="0.7"/><circle cx="195" cy="90" r="5" fill="#ef4444" opacity="0.7"/>
  <circle cx="220" cy="210" r="5" fill="#ef4444" opacity="0.7"/><circle cx="240" cy="130" r="5" fill="#ef4444" opacity="0.7"/>
  <circle cx="260" cy="160" r="5" fill="#ef4444" opacity="0.7"/><circle cx="280" cy="95" r="5" fill="#ef4444" opacity="0.7"/>
  <!-- Spurious trend line -->
  <line x1="70" y1="210" x2="285" y2="100" stroke="#ef4444" stroke-width="2" stroke-dasharray="6,4" opacity="0.6"/>
  <text x="175" y="288" text-anchor="middle" font-size="10" fill="#dc2626" font-style="italic">Apparent trend is just random scatter</text>
  <!-- Right panel: clear signal with many points -->
  <rect x="360" y="42" width="290" height="240" rx="6" fill="#f8fafc" stroke="#e2e8f0"/>
  <text x="505" y="62" text-anchor="middle" font-size="12" font-weight="600" fill="#16a34a">Larger Sample (n=100) — Pattern Is Clearer</text>
  <!-- Axes -->
  <line x1="390" y1="250" x2="625" y2="250" stroke="#9ca3af" stroke-width="1.5"/>
  <line x1="390" y1="250" x2="390" y2="75" stroke="#9ca3af" stroke-width="1.5"/>
  <text x="505" y="272" text-anchor="middle" font-size="11" fill="#6b7280">Measurement</text>
  <!-- Many points clustered around trend -->
  <circle cx="400" cy="225" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="405" cy="218" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="410" cy="230" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="415" cy="210" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="420" cy="222" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="425" cy="215" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="430" cy="205" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="435" cy="212" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="440" cy="200" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="445" cy="208" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="450" cy="195" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="455" cy="190" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="460" cy="198" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="465" cy="185" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="470" cy="192" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="475" cy="180" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="480" cy="188" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="485" cy="175" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="490" cy="182" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="495" cy="170" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="500" cy="178" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="505" cy="165" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="510" cy="172" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="515" cy="160" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="520" cy="168" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="525" cy="155" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="530" cy="162" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="535" cy="150" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="540" cy="158" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="545" cy="145" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="550" cy="152" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="555" cy="140" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="560" cy="148" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="565" cy="135" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="570" cy="142" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="575" cy="128" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="580" cy="138" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="585" cy="122" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="590" cy="130" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="595" cy="118" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="600" cy="125" r="3" fill="#3b82f6" opacity="0.4"/><circle cx="605" cy="110" r="3" fill="#3b82f6" opacity="0.4"/>
  <circle cx="610" cy="115" r="3" fill="#3b82f6" opacity="0.4"/>
  <!-- Clear trend line -->
  <line x1="395" y1="232" x2="615" y2="108" stroke="#2563eb" stroke-width="2.5" opacity="0.85"/>
  <text x="505" y="288" text-anchor="middle" font-size="10" fill="#16a34a" font-style="italic">More data reduces sampling noise; it does not remove bias</text>
</svg>
</div>

Consider a simple experiment. You flip a coin 10 times and get 8 heads. Is the coin biased? Your intuition says maybe — 8 out of 10 is a lot of heads. But if you do the math, a fair coin produces 8 or more heads about 5.5% of the time. That is unlikely, but not astronomically so. You might just be unlucky.

Now flip the coin 100 times and get 80 heads. Under a fair-coin model, the probability of 80 or more heads is about 0.00000000056—roughly one in 1.8 billion. That is strong evidence against the fair-coin model, although it does not diagnose why the result occurred.

The observed proportion (80% heads) is the same in both cases. What changed is the **sample size**. Under the fair-coin model, 8 of 10 can occur occasionally; 80 of 100 is much harder to explain by sampling variation alone.

The same reasoning applies to the illustrative cancer trial: a small early study leaves more room for sampling variation. A larger, well-designed study usually narrows that uncertainty, but size alone does not remove bias, measurement problems, or differences between study populations.

> **Common pitfall:** Small studies frequently produce dramatic-looking results. This is not because small studies discover larger effects — it is because small samples are inherently noisy, and noise occasionally looks like a big signal. This phenomenon is called the "winner's curse" and it haunts biomedical research.

## The Cost of Being Wrong

In a binary hypothesis-testing decision, two error types are possible, and they have different consequences.

### Type I Error: The False Alarm

A Type I error occurs when you conclude there is an effect when there is none. You declare the coin biased when it is actually fair. You approve a drug that does not work.

A laboratory example is advancing a compound because one noisy screen crosses a threshold even though the compound has no repeatable effect. Replication, quality control, and multiplicity correction reduce this risk; no single threshold eliminates it.

### Type II Error: The Missed Discovery

A Type II error occurs when the test fails to reject its null model even though the specified alternative is true. You miss evidence that the coin is biased or fail to detect a treatment effect the study was designed to find.

A laboratory example is discarding a promising treatment because a small, noisy experiment could not distinguish its effect from background variation. A non-significant result may reflect no important effect, inadequate precision, poor measurement, or a mismatched analysis; the interval and design help distinguish these possibilities.

The difficulty is that a missed effect is not visible in the final table. A wide confidence interval may include both no meaningful benefit and a worthwhile benefit. That is why planning should consider the smallest effect that matters, expected variation, attrition, and the experimental unit before data collection. Day 18 develops power and sample-size planning.

| Error Type | What Happens | Consequence | Biology Example |
|---|---|---|---|
| Type I (False Positive) | Conclude effect exists when it does not | Wasted resources, patient harm | Approving ineffective drug |
| Type II (False Negative) | Miss a real effect | Lost discoveries, delayed treatments | Rejecting H. pylori hypothesis |
| Correct non-rejection | Test does not reject when H0 is true | Avoid an unnecessary claim | Do not advance an inactive screening hit |
| Correct detection | Reject H0 when the specified alternative is true | Detect a planned signal | Detect a repeatable treatment effect |

<div style="text-align: center; margin: 2em 0;">
<svg width="660" height="310" viewBox="0 0 660 310" xmlns="http://www.w3.org/2000/svg" style="background: #fafbfc; border: 1px solid #e5e7eb; border-radius: 8px;">
  <text x="330" y="28" text-anchor="middle" font-size="15" font-weight="bold" fill="#1e293b">Decision Outcomes: The 2x2 Reality</text>
  <!-- Column headers -->
  <text x="370" y="60" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">Reality: No Effect (H0 true)</text>
  <text x="560" y="60" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">Reality: Effect Exists (H1 true)</text>
  <!-- Row headers -->
  <text x="120" y="132" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">Decision:</text>
  <text x="120" y="148" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">"Significant"</text>
  <text x="120" y="232" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">Decision:</text>
  <text x="120" y="248" text-anchor="middle" font-size="12" font-weight="600" fill="#6b7280">"Not Significant"</text>
  <!-- Grid lines -->
  <line x1="200" y1="70" x2="200" y2="290" stroke="#d1d5db" stroke-width="1.5"/>
  <line x1="470" y1="70" x2="470" y2="290" stroke="#d1d5db" stroke-width="1.5"/>
  <line x1="200" y1="70" x2="640" y2="70" stroke="#d1d5db" stroke-width="1.5"/>
  <line x1="200" y1="180" x2="640" y2="180" stroke="#d1d5db" stroke-width="1.5"/>
  <line x1="200" y1="290" x2="640" y2="290" stroke="#d1d5db" stroke-width="1.5"/>
  <line x1="640" y1="70" x2="640" y2="290" stroke="#d1d5db" stroke-width="1.5"/>
  <!-- Type I Error (False Positive) -->
  <rect x="202" y="72" width="266" height="106" fill="#fef2f2" rx="4"/>
  <text x="335" y="110" text-anchor="middle" font-size="14" font-weight="bold" fill="#dc2626">TYPE I ERROR</text>
  <text x="335" y="130" text-anchor="middle" font-size="12" fill="#dc2626">False Positive (alpha)</text>
  <text x="335" y="150" text-anchor="middle" font-size="11" fill="#7f1d1d">Approve drug that doesn't work</text>
  <text x="335" y="166" text-anchor="middle" font-size="11" fill="#7f1d1d">Chosen alpha controls this rate</text>
  <!-- True Positive -->
  <rect x="472" y="72" width="166" height="106" fill="#f0fdf4" rx="4"/>
  <text x="555" y="110" text-anchor="middle" font-size="14" font-weight="bold" fill="#16a34a">CORRECT</text>
  <text x="555" y="130" text-anchor="middle" font-size="12" fill="#16a34a">True Positive (Power)</text>
  <text x="555" y="150" text-anchor="middle" font-size="11" fill="#14532d">Detect real treatment</text>
  <text x="555" y="166" text-anchor="middle" font-size="11" fill="#14532d">Power depends on effect and design</text>
  <!-- True Negative -->
  <rect x="202" y="182" width="266" height="106" fill="#f0fdf4" rx="4"/>
  <text x="335" y="220" text-anchor="middle" font-size="14" font-weight="bold" fill="#16a34a">CORRECT</text>
  <text x="335" y="240" text-anchor="middle" font-size="12" fill="#16a34a">True Negative</text>
  <text x="335" y="260" text-anchor="middle" font-size="11" fill="#14532d">Do not claim an unsupported effect</text>
  <text x="335" y="276" text-anchor="middle" font-size="11" fill="#14532d">Probability = 1 - alpha</text>
  <!-- Type II Error (False Negative) -->
  <rect x="472" y="182" width="166" height="106" fill="#fef2f2" rx="4"/>
  <text x="555" y="220" text-anchor="middle" font-size="14" font-weight="bold" fill="#dc2626">TYPE II ERROR</text>
  <text x="555" y="240" text-anchor="middle" font-size="12" fill="#dc2626">False Negative (beta)</text>
  <text x="555" y="260" text-anchor="middle" font-size="11" fill="#7f1d1d">Miss a real treatment</text>
  <text x="555" y="276" text-anchor="middle" font-size="11" fill="#7f1d1d">Rate = 1 - Power</text>
</svg>
</div>

> **Clinical relevance:** In diagnostic testing, Type I errors produce false positives (telling a healthy person they have cancer) and Type II errors produce false negatives (telling a cancer patient they are healthy). Both are harmful, but in different ways. The balance between them is one of the central tensions in medicine.

## Why Biological Data Is Challenging

Biological measurements combine real variation among organisms with technical variation from sampling and instruments. Three features appear repeatedly.

### 1. Biological Variability

Organisms, tissues, cells, and molecules vary. Even genetically similar organisms in controlled conditions can differ in expression, immune response, growth, and other outcomes. This is biological variation, not automatically experimental error.

One observation may be important, but it cannot by itself describe repeatability or population variation. Independent biological replication helps separate an individual observation from a reproducible pattern.

### 2. Measurement Noise

RNA-seq quantification depends on sample handling, library preparation, read depth, alignment, and normalization. Mass-spectrometry results depend on preparation, calibration, detection, and ionization. Some errors are random; others are systematic.

Statistics can quantify variation and compare explanations, but it cannot always identify “true signal” from values alone. Controls, calibration, replication, and recorded metadata are essential.

### 3. Massive Parallel Testing

Modern biology is high-dimensional. A genome-wide association study (GWAS) tests millions of genetic variants. A differential expression analysis tests 20,000 genes. A proteomics experiment quantifies thousands of proteins. A drug screen tests hundreds of compounds.

If 20,000 valid null hypotheses were all true and each were tested at 0.05, the expected number crossing that unadjusted threshold would be about 1,000. This is why genomics workflows must account for multiplicity.

To make this concrete, consider a differential expression analysis. You measure expression of 20,000 genes in treatment versus control. Even if the treatment does absolutely nothing — affects zero genes — testing each gene at &alpha; = 0.05 will flag approximately 1,000 genes as "significant." If you published a paper claiming these 1,000 genes are treatment-responsive, every single one would be a false positive.

Multiple-testing methods control a chosen error rate across the collection of questions. Day 12 explains family-wise error and false-discovery-rate approaches; the appropriate choice depends on whether the goal is strict confirmation or discovery followed by validation.

> **Key insight:** Biological variation, technical effects, and thousands of simultaneous measurements make design, quality control, uncertainty, and multiplicity part of the science—not administrative details.

## A Tour of What Lies Ahead

This book provides a working foundation for analysing biological data and communicating uncertainty. It does not replace specialist training or collaboration for high-stakes studies. Here is the learning path:

<div style="text-align: center; margin: 2em 0;">
<svg width="680" height="200" viewBox="0 0 680 200" xmlns="http://www.w3.org/2000/svg" style="background: #fafbfc; border: 1px solid #e5e7eb; border-radius: 8px;">
  <text x="340" y="24" text-anchor="middle" font-size="15" font-weight="bold" fill="#1e293b">30-Day Roadmap</text>
  <!-- Timeline line -->
  <line x1="40" y1="80" x2="640" y2="80" stroke="#d1d5db" stroke-width="3" stroke-linecap="round"/>
  <!-- Week 1 -->
  <rect x="40" y="65" width="100" height="30" rx="15" fill="#2563eb"/>
  <text x="90" y="85" text-anchor="middle" font-size="11" font-weight="600" fill="white">Week 1</text>
  <text x="90" y="112" text-anchor="middle" font-size="10" fill="#2563eb" font-weight="600">Foundations</text>
  <text x="90" y="126" text-anchor="middle" font-size="9" fill="#6b7280">Descriptive stats,</text>
  <text x="90" y="138" text-anchor="middle" font-size="9" fill="#6b7280">distributions, prob.</text>
  <!-- Week 2 -->
  <rect x="160" y="65" width="100" height="30" rx="15" fill="#7c3aed"/>
  <text x="210" y="85" text-anchor="middle" font-size="11" font-weight="600" fill="white">Week 2</text>
  <text x="210" y="112" text-anchor="middle" font-size="10" fill="#7c3aed" font-weight="600">Hypothesis Testing</text>
  <text x="210" y="126" text-anchor="middle" font-size="9" fill="#6b7280">CIs, p-values, t-tests,</text>
  <text x="210" y="138" text-anchor="middle" font-size="9" fill="#6b7280">non-parametric</text>
  <!-- Week 3 -->
  <rect x="280" y="65" width="100" height="30" rx="15" fill="#16a34a"/>
  <text x="330" y="85" text-anchor="middle" font-size="11" font-weight="600" fill="white">Week 3</text>
  <text x="330" y="112" text-anchor="middle" font-size="10" fill="#16a34a" font-weight="600">Beyond Two Groups</text>
  <text x="330" y="126" text-anchor="middle" font-size="9" fill="#6b7280">ANOVA, chi-square,</text>
  <text x="330" y="138" text-anchor="middle" font-size="9" fill="#6b7280">regression, MTC</text>
  <!-- Week 4 -->
  <rect x="400" y="65" width="100" height="30" rx="15" fill="#dc2626"/>
  <text x="450" y="85" text-anchor="middle" font-size="11" font-weight="600" fill="white">Week 4</text>
  <text x="450" y="112" text-anchor="middle" font-size="10" fill="#dc2626" font-weight="600">Advanced Methods</text>
  <text x="450" y="126" text-anchor="middle" font-size="9" fill="#6b7280">Survival, logistic reg,</text>
  <text x="450" y="138" text-anchor="middle" font-size="9" fill="#6b7280">PCA, clustering</text>
  <!-- Week 5-6 -->
  <rect x="520" y="65" width="120" height="30" rx="15" fill="#ea580c"/>
  <text x="580" y="85" text-anchor="middle" font-size="11" font-weight="600" fill="white">Weeks 5-6</text>
  <text x="580" y="112" text-anchor="middle" font-size="10" fill="#ea580c" font-weight="600">Genomics & Practice</text>
  <text x="580" y="126" text-anchor="middle" font-size="9" fill="#6b7280">DE analysis, ML,</text>
  <text x="580" y="138" text-anchor="middle" font-size="9" fill="#6b7280">study design, capstone</text>
  <!-- Milestone markers -->
  <circle cx="210" cy="160" r="4" fill="#7c3aed"/>
  <text x="210" y="175" text-anchor="middle" font-size="9" fill="#7c3aed">Day 8: Compare groups</text>
  <circle cx="450" cy="160" r="4" fill="#dc2626"/>
  <text x="450" y="175" text-anchor="middle" font-size="9" fill="#dc2626">Day 17: Survival curves</text>
  <circle cx="580" cy="160" r="4" fill="#ea580c"/>
  <text x="580" y="175" text-anchor="middle" font-size="9" fill="#ea580c">Day 30: Full GWAS plan</text>
</svg>
</div>

**Week 1 (Days 1-5): Foundations.** You will learn to summarize data, understand distributions, reason about probability, and appreciate why sample size matters. These are the tools you need before you can test any hypothesis.

**Week 2 (Days 6-10): Hypothesis Testing.** You will learn confidence intervals, p-values, t-tests, and non-parametric alternatives. By Day 8, you will be able to rigorously determine whether two groups differ. By Day 10, you will know when to use (and when to avoid) parametric tests.

**Week 3 (Days 11-15): Beyond Two Groups.** ANOVA, chi-square tests, correlation, regression, and multiple testing correction. You will analyze multi-group experiments, categorical data, and learn why "correlation does not imply causation" is more nuanced than it sounds.

**Week 4 (Days 16-20): Advanced Methods.** Survival analysis, logistic regression, principal component analysis, and clustering. You will build Kaplan-Meier curves, classify patients, reduce high-dimensional data, and find natural groupings in gene expression datasets.

**Week 5 (Days 21-25): Genomics Applications.** Differential expression analysis, enrichment testing, multiple testing correction in practice, and Bayesian thinking. The methods that power modern computational biology.

**Week 6 (Days 26-30): Real-World Practice.** Power analysis and study design, meta-analysis, machine learning basics, reproducible research practices, and a capstone project that ties everything together.

By Day 8, you will know if two groups truly differ. By Day 17, you will build survival curves that predict patient outcomes. By Day 22, you will analyze differential gene expression. By Day 30, you will design a complete statistical analysis plan for a GWAS.

Each day follows the same pattern: a real-world problem that motivates the method, the conceptual framework, hands-on BioLang code, comparisons with Python and R, and exercises to cement your understanding. The emphasis throughout is on **understanding** — not memorizing formulas, but developing the intuition to know which method to use and why.

## The Statistician's Mindset

Before we dive into formulas and code, internalize these four questions. Ask them every time you look at data, read a paper, or plan an experiment:

### 1. How variable is it?

A mean without a measure of spread is almost meaningless. "Average tumor size decreased by 2 cm" sounds impressive until you learn that the standard deviation was 4 cm. Always ask: what is the spread?

### 2. Could chance explain this?

The human brain is wired to see patterns, even in random noise. We see faces in clouds, constellations in random stars, and trends in stock prices. Before accepting any pattern as real, quantify the probability that it arose by chance. This is the essence of hypothesis testing.

### 3. How big is the effect?

Statistical significance and practical significance are not the same thing. With a large enough sample, you can detect a difference of 0.001 grams in tumor weight with p < 0.001. But is a one-milligram difference clinically meaningful? Always report effect sizes alongside p-values.

### 4. Is my sample representative?

If you study the genetics of heart disease using only patients from a single hospital in Boston, your results may not generalize to patients in rural India. If you select only the "best" cell lines for your experiment, your conclusions may not extend to primary cells. Sampling bias is the silent assassin of biomedical research.

### Putting the Mindset into Practice

These four questions are not abstract philosophy. They are a practical checklist:

| Question | When Reading a Paper | When Designing an Experiment |
|---|---|---|
| How variable? | Check SD, IQR, range | Plan enough replicates |
| Could chance explain it? | Scrutinize p-values, CI | Pre-register analysis plan |
| How big is the effect? | Look for effect sizes, not just significance | Define minimum meaningful difference |
| Representative sample? | Check inclusion criteria, demographics | Match your sample to target population |

You will encounter these questions again and again throughout this book. By Day 30, they will be second nature — the automatic mental checklist of a statistically literate scientist.

> **Key insight:** Statistics is not a set of tests to run after the experiment. It is a way of thinking that should inform every stage — from study design to data collection to analysis to interpretation. The best time to consult a statistician is before you collect a single data point.

## The Burden of Proof

In everyday life, we make decisions based on intuition, anecdote, and authority. "My grandmother smoked until 95, so smoking cannot be that bad." "This supplement worked for my friend, so it must be effective." "The famous professor says this treatment works, so it must."

Science demands a higher standard. The burden of proof rests on the claimant. If you claim a drug works, you must demonstrate it with evidence strong enough to withstand scrutiny. If you claim a gene is associated with a disease, you must show that the association is unlikely to be a coincidence.

Statistics provides the machinery for this burden of proof. It forces you to be explicit about your assumptions, quantify your uncertainty, and acknowledge the limits of your data. It is, in essence, formalized humility.

Consider the claim "Vitamin D supplements reduce cancer risk." An anecdote is worthless: your uncle took vitamin D and did not get cancer. A small observational study is weak: 50 people who took vitamin D had fewer cancers than 50 who did not — but maybe the vitamin D group was healthier to begin with (confounding). A large randomized controlled trial with 25,000 participants, pre-registered outcomes, and proper statistical analysis is strong evidence. Each step up the ladder requires more statistical sophistication.

The hierarchy of evidence is, fundamentally, a hierarchy of statistical rigor:

| Evidence Level | Design | Statistical Rigor |
|---|---|---|
| Weakest | Case report / anecdote | None |
| Weak | Case series | Descriptive only |
| Moderate | Observational study | Potential confounding |
| Strong | Randomized controlled trial | Causal inference possible |
| Strongest | Meta-analysis of multiple RCTs | Pooled estimates, high power |

This book will equip you to evaluate and produce evidence at every level of this hierarchy.

## The Numbers Tell a Story

To bring this all together, let us look at a real-world scenario that illustrates every concept from today.

A research group publishes a paper claiming that a new biomarker predicts response to immunotherapy. Their study: 24 patients (12 responders, 12 non-responders). They measure the biomarker level in each patient and find a "statistically significant" difference (p = 0.03).

Here is what a statistical thinker would ask:

**How variable is it?** The biomarker levels range from 2 to 200 ng/mL. The standard deviation within each group is enormous — nearly as large as the difference between groups. The signal is weak relative to the noise.

**Could chance explain it?** With only 12 per group and high variability, the p-value of 0.03 is fragile. If you removed two extreme patients, it becomes 0.12. The result is not robust.

**How big is the effect?** The difference in medians is 15 ng/mL, but the overlap between groups is substantial. Many responders have lower biomarker levels than many non-responders. The effect size (Cohen's d) is only 0.4 — a "small to medium" effect.

**Is the sample representative?** All patients came from a single institution, were predominantly male, and had a specific tumor subtype. Whether the biomarker works in a broader population is unknown.

A naive reader sees "p < 0.05, significant." A statistically literate reader sees a fragile, underpowered result from a non-representative sample with a modest effect size. These are different conclusions from the same data.

## A Preview of the Tools

Throughout this book, you will use BioLang to perform statistical analyses. Here is a tiny glimpse of what Day 2 will look like — just to whet your appetite:

```
# Tomorrow, you'll summarize 10,000 quality scores in one line:
# let stats = summary(quality_scores)
#
# And visualize them instantly:
# histogram(quality_scores, {bins: 50, title: "Sequencing Quality Distribution"})
```

But today is about the **why**, not the **how**. The tools are only as good as the thinking behind them. A researcher who understands why a t-test exists will use it correctly even with imperfect software. A researcher who merely knows how to call a t-test function will misuse it regularly, regardless of how elegant the software is.

## Exercises

### Exercise 1: The Newspaper Test

Find a news article reporting a scientific or medical finding (e.g., "Coffee reduces cancer risk by 15%"). Write down your answers to these four questions:

- (a) What was the sample size? If the article does not mention it, what does that tell you?
- (b) Could the result be due to chance? What would you need to know to answer this?
- (c) Is the effect size meaningful in practice? A 2% reduction in cancer risk sounds different from a 50% reduction.
- (d) Is the sample representative of the population you care about? Who was studied, and who was not?

If the article does not provide enough information to answer these questions, that itself is informative. Most science journalism omits sample sizes, effect sizes, and confidence intervals — precisely the information you need to evaluate the claim.

### Exercise 2: Coin Flip Thought Experiment

Without doing any math, estimate the following:

- If you flip a fair coin 20 times, what is the probability of getting exactly 10 heads?
- What about 15 or more heads out of 20?
- What about 20 heads in a row?

Write down your guesses. We will revisit this on Day 4 with the tools to compute exact answers, and you can see how well your intuition calibrated.

### Exercise 3: Reproducibility Reflection

Think about a result from your own work (or a paper you have read) that you found surprising or striking. List three reasons why the result might fail to reproduce if someone repeated the experiment. For each reason, identify whether it is:

- Statistical (sample size, random variation, multiple testing)
- Methodological (different protocols, reagent lots, equipment)
- Biological (different cell lines, patient populations, environmental conditions)

### Exercise 4: Type I vs Type II in Your Field

Identify one example each of a Type I error (false positive) and a Type II error (false negative) that would be particularly damaging in your area of biology. For each:

- Describe the scenario concretely
- Estimate the consequences (financial, clinical, scientific)
- State which error type you consider more dangerous in your context, and why

### Exercise 5: Spotting P-Hacking

A paper reports testing a drug on patients across 8 different cancer subtypes. Only one subtype shows a significant result (p = 0.04). The paper's title highlights this positive finding. What statistical concerns should this raise? How many false positives would you expect by chance when testing 8 subtypes at &alpha; = 0.05?

## Key Takeaways

- Statistics is the science of learning from data in the presence of uncertainty — it is essential, not optional, for biological research.
- Results may fail to reproduce because of sampling variation, bias, protocol differences, selective reporting, confounding, or analysis choices; transparent workflows make these causes easier to examine.
- Signal vs. noise is the fundamental statistical question: the same percentage difference can be meaningful or meaningless depending on sample size.
- Type I errors (false positives) waste resources and can cause harm; Type II errors (false negatives) cause missed discoveries. Neither can be eliminated — only managed.
- Biological analysis must handle biological variation, technical effects, and often thousands of simultaneous measurements.
- The statistician's mindset asks four questions: How variable? Could chance explain it? How big is the effect? Is the sample representative?
- Statistics should inform every stage of research, from design through interpretation — not just the analysis phase.

## What's Next

Tomorrow, we roll up our sleeves and meet data. You will summarize 10,000 numbers with means, medians, standard deviations, quantiles, and visual checks. You will see how extreme observations affect the mean, how box plots and histograms show complementary evidence, and how BioLang can organize QC clues without making the scientific accept/reject decision for you. Day 2 is where the hands-on work begins.
