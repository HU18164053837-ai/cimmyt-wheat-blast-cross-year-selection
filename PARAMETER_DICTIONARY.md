# Locked parameter dictionary

## Identity and analysis units

- Cross-year material identity uses a nonmissing CIMMYT GID only.
- A nursery-environment stratum is nursery x location x year x sowing.
- Unobserved genotype-environment combinations are structural absences and are not imputed.
- The exact 15HZAN duplicate of 14HZAN is excluded from the primary long table.

## Percentile definitions

- Material percentile within a nursery-environment stratum is `rank(disease, ties.method="average") / n`. Lower values are favorable.
- Environment-component percentile rank is `(rank(x, ties.method="average") - 1) / (N - 1)`. It maps the lowest component value to 0 and the highest to 1.
- These two mappings serve different estimands and are not interchangeable.

## Cross-year eligibility

All conditions are required.

- Nonmissing GID.
- At least 2 observed environmental years.
- At least 2 observed nurseries.
- At least 4 phenotype observations.

## Robustness tiers

Strong requires all of the following.

- Median material percentile <= 0.40.
- 75th percentile of material percentile <= 0.45.
- Worst-year median material percentile <= 0.50.
- 90th percentile disease <= 10%.
- Maximum observed disease <= 20%.

Moderate requires all of the following and is assigned only when the strong rule is not met.

- Median material percentile <= 0.45.
- 75th percentile of material percentile <= 0.60.
- Worst-year median material percentile <= 0.65.
- 90th percentile disease <= 20%.
- Maximum observed disease <= 50%.

All other eligible GIDs are classified as not robust.

## Priority evidence set

- Two-year evidence enters the priority set only when the tier is strong.
- Three- to four-year evidence enters when the tier is strong or moderate.
- Five- to seven-year evidence enters when the tier is strong or moderate.
- The union of these deterministic conditions produces the 15-GID priority evidence set.

## Pedigree de-duplication

- Pedigree strings are tokenized into named components.
- Pedigree strings are converted to token sets, pairwise Jaccard distances are calculated as `1 - similarity`, and average-linkage hierarchical clustering is cut at distance 0.70 (corresponding to a nominal similarity level of 0.30) to define the primary exploratory pedigree clusters.
- Long-term and intermediate candidates are retained as evidence anchors.
- Within a redundant two-year strong cluster, candidates are ordered by lower disease 90th percentile, lower maximum disease, and lower median material percentile, in that order.
- One representative is retained per cluster. The exact-pedigree pair 7627645/7627655 retains GID 7627645.
- This step produces the 14-GID nonredundant resistance portfolio.

These thresholds are operational reference rules rather than preregistered biological cutoffs. They jointly require relative-rank evidence and control of the raw disease tail. The accompanying threshold grid quantifies how the candidate count changes around the reference settings; it does not establish that the resulting 14-member set is uniquely optimal.

## Rolling temporal evaluation

- For cutoff year t, training uses observations through t and evaluation uses t+1.
- The evaluation population is restricted to GIDs reobserved in t+1 with at least 3 evaluation observations.
- Training eligibility requires at least 3 training observations and at least 2 training locations.
- Training selection requires median percentile <= 0.40, 75th percentile <= 0.50, disease 90th percentile <= 20%, and a negative LOCAL CHECK contrast when available.
- A conditional next-year pass requires at least 3 observations, median percentile <= 0.50, 75th percentile <= 0.65, disease 90th percentile <= 30%, and no unfavorable reversal against the check when available.
- The analysis evaluates a simplified historical signal rule among reobserved GIDs. It does not validate the complete 14-member portfolio construction algorithm.

## Environment information score

- Pressure, discrimination, nursery coverage, and nonsaturation each receive weight 0.25.
- Nonsaturation is `1 - percentile_rank(floor_percentage + ceiling_percentage)`.
- The score is a historical archive information profile for prospective trial design, not a causal environmental effect or independently validated site-quality metric.

## Uncertainty and interpretation

- LOCAL CHECK intervals use 2,000 complete environment-cluster bootstrap samples.
- Background-adjusted intervals use 1,000 complete environment-cluster bootstrap samples.
- Temporal uncertainty uses 10,000 transition-cluster and GID-cluster bootstrap samples and 10,000 within-transition permutations with plus-one correction.
- Candidate intervals are selection-conditional descriptive intervals. They are not independent confirmation of resistance.

## Full-funnel analyst-choice sensitivity

- The 768 strong-rule combinations are paired with the same strong-to-moderate relaxation used by the reference rule.
- Each paired rule is passed through the locked evidence-duration gates.
- Pedigree de-duplication is repeated at nominal Jaccard-similarity levels 0.20, 0.30, and 0.40, corresponding to average-linkage distance cuts 0.80, 0.70, and 0.60.
- The resulting 2,304 complete decision paths produce a final-membership frequency for every eligible GID.
- Frequency >=0.80 is labeled stable core, 0.50-0.799 moderately stable, and <0.50 choice-sensitive. These labels describe analyst-choice sensitivity and are not biological probabilities.

## Analysis-unit hierarchy

- 38 location-year-sowing environments are used for equal-nursery environmental summaries.
- 148 nursery-environment strata preserve nursery-specific composition within those broader environments.
- A candidate by nursery-environment event is the within-cell median for one focal GID in one nursery-environment stratum. The background-adjusted analysis contains 275 observed events.
- Background-adjusted intervals resample the complete nursery-environment cluster key, defined as nursery plus location-year-sowing environment.
