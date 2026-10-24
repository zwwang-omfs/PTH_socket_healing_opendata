# Two-way ANOVA of in vitro endpoints

## Scope

- Input: `data/statistics/in_vitro`
- Endpoints: 12
- Analysis unit: each CSV row is one independent experimental observation.
- No technical-replicate averaging or batch pairing was applied.
- The source CSV files were not modified.

## Model

For every endpoint:

`value ~ categorical 4-HC concentration * teriparatide`

Type III tests use sum-to-zero contrasts. Brown-Forsythe testing across
factorial cells was used as a prespecified covariance rule: classical
covariance when p >= 0.05 and HC3 heteroskedasticity-robust covariance when
p < 0.05. HC3 was selected for 0 of 12 endpoints.

## Planned comparisons

1. Teriparatide versus no teriparatide at 0 uM: no multiplicity correction.
2. Each nonzero 4-HC concentration versus 0 uM without teriparatide:
   Dunnett single-step multivariate-t correction.
3. Teriparatide versus no teriparatide at each nonzero concentration:
   one rescue family with single-step multivariate-t correction.
4. Each combined-treatment cell versus 0 uM/no-teriparatide control:
   a separate exploratory family with single-step multivariate-t correction.

The multivariate-t calculations use reproducible scrambled Sobol
quasi-Monte Carlo integration (`seed=20260615`, `2^19` draws per
family). Simultaneous 95% confidence intervals are reported for adjusted
families.

## Rescue interpretation

- `partial rescue`: adjusted rescue comparison is significant in the
  direction of baseline and the adjusted diagonal comparison still differs
  from baseline in the original impairment direction.
- `improvement; complete rescue not established`: rescue is significant
  toward baseline, while the diagonal comparison is not significant.
- `complete rescue` is never concluded here because no equivalence margin
  was prespecified. A nonsignificant difference from baseline does not
  establish equivalence.

## Source column labels

Some source files use the outcome column label `OD450` despite filenames
indicating scratch or tube-network metrics. Results retain both the filename
endpoint and original source column so these labels can be checked before
manuscript use.
