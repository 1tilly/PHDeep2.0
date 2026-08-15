# Stage Contracts

## Purpose
This document defines the stable input/output contracts for the core PHDeep2.0 stages.
These contracts are the compatibility surface for local and batch backends.

## Stage: `bed_to_training`
Description:
Transforms ENCODE-style metadata + BED inputs into a feature-labeled training regions file.

Inputs:
- `metadata.tsv` with required columns:
  - `File accession`
  - `File format`
  - `File type`
  - `Experiment accession`
  - `Output type`
  - `Assay`
  - `Biosample term name`
  - `Experiment target`
  - `File download URL`
  - `Biological replicate(s)`
  - `Technical replicate(s)`
  - `Biosample treatments`
  - `File assembly`
  - `File analysis title`
  - `Biosample organism`
- `input_bed_files/` containing `{File accession}.bed.gz`

Outputs:
- `training_regions.bed` columns:
  - `chrom`
  - `start`
  - `end`
  - `feature`
- `features.txt`: one feature per line.
- `read_errors.txt`: parser/read errors.

Guarantees:
- Output `chrom` values are normalized to `chr*`.
- Output rows are limited to chromosomes `chr1-22`, `chrX`, `chrY`.

## Stage: `train` (canonical first path)
Description:
Trains one model on fixed-length sequence windows from training regions.

Inputs:
- `training_regions.bed`
- reference FASTA
- model/training config

Outputs:
- model checkpoint artifact
- run metrics (at minimum: loss, AUROC where available)
- run metadata (config + seed)

## Stage: `predict`
Modes:
- `reference`: predict scores for reference windows.
- `variant`: predict ref/alt effects for variant records.

Inputs:
- model checkpoint
- one of:
  - reference regions bed (`reference` mode)
  - variants feather/tsv (`variant` mode)

Outputs (`mode=reference`, via `ReferencePredictor.predict_bed`):
- table with columns:
  - `chrom`, `start`, `end`
  - `pred_<i>` for each model output feature `i`

Outputs (`mode=variant`, via `VariantEffectPredictor.predict_variants`):
- the input variant columns (`chromosome`, `start`, `reference`, `alternate`,
  plus any caller-supplied metadata columns such as a grouping column) with
  three columns appended per model output feature `i`:
  - `ref_pred_<i>`, `alt_pred_<i>`, `delta_<i>` (= `alt_pred_<i> - ref_pred_<i>`)
- `config.pipeline_config.PREDICTION_OUTPUT_COLUMNS` captures only the
  fixed `(chromosome, start, reference, alternate)` columns of the
  `variant` mode output — a fixed-width tuple can't express "three columns
  per model feature", so the `ref_pred_*`/`alt_pred_*`/`delta_*` columns
  are not enumerated there.

## Stage: `aggregate`
Description:
Turns per-variant `predict` (mode=variant) output into the per-variant
weights table SKAT-O needs, via `build_variant_weights_table`. This is
**not** collapsed to one row per group/gene — it is one row per variant,
canonically sorted by `(chromosome, start, end, alternate)`, with a
`group` column identifying which gene/region each variant belongs to.

Inputs:
- one or more `predict` (mode=variant) output tables, all with matching
  `delta_*`/`ref_pred_*`/`alt_pred_*` columns (the runner rejects mismatched
  schemas across input files rather than silently NaN-filling them via
  `pd.concat`)
- `group_col`: which input column to use for grouping (e.g. `gene_symbol`)

Outputs:
- `output_scores` (required): the per-variant weights table.
  `config.pipeline_config.AGGREGATION_OUTPUT_COLUMNS` is the exact,
  ordered column list (verified against `build_variant_weights_table`'s
  real emission order by the PH2-013 golden fixture, not just a set):
  `chromosome, start, end, reference, alternate, variant_id, group,
  n_features, eis_ref, eis_alt, eis_diff, abs_delta_max, abs_delta_sum,
  l2_delta`. `eis_ref`/`eis_alt` are per-variant sums of ref/alt
  predictions across model features; `eis_diff` is their difference;
  `abs_delta_max`/`abs_delta_sum`/`l2_delta` summarize the per-feature
  delta vector. Rows with a null `group_col` value are retained (not
  dropped), with `group` set to null.
- `output_genotypes` (optional, requires `sample_ids_file`): a variant x
  sample genotype (dosage) matrix, via `build_genotype_matrix`, in the
  same variant order as `output_scores`. Columns:
  `config.pipeline_config.GENOTYPE_MATRIX_KEY_COLUMN` (`variant_id`) plus
  one integer column per sample id, values in `{0, 1, 2, 9}` (9 = missing
  genotype call).
- `output_group_summary` (optional): the older one-row-per-group summary
  from `aggregate_variant_scores`, kept for callers that still want a
  group-level rollup rather than per-variant weights.

## Stage: `stats`
Description:
Runs SKAT-O per group (gene/region) on a real genotype matrix and a real
per-sample phenotype, via `run_skat_o_from_feather`.

Inputs:
- `input_scores`: the `aggregate` stage's `output_scores` (per-variant
  weights table)
- `input_genotypes`: the `aggregate` stage's `output_genotypes` (variant x
  sample dosage matrix)
- `phenotype_table` (**required**, cohort input this pipeline does not
  generate): a feather or TSV file with a sample-id column
  (`sample_id_col`, default `sample_id`) and a phenotype column
  (`phenotype_col`, default `phenotype`, values `1`=case/`0`=control).
  Real per-sample phenotypes (and, for a full analysis, covariates and
  kinship) must come from the cohort's own data management — nothing in
  this pipeline synthesizes them.
- `group_col` (default `group`), `weight_col` (optional per-variant weight
  column from `input_scores`, e.g. `eis_diff` or `abs_delta_max`),
  `min_variants` (groups with fewer variants are skipped), `method`
  (default `"optimal.adj"`)

Outputs:
- result table with columns (`config.pipeline_config.STATS_OUTPUT_COLUMNS`,
  an exact ordered list — verified against `run_skat_o`'s real emission
  order by the PH2-013 golden fixture, not just a set): `feature_id,
  n_variants, n_samples, p_value, p_value_burden, p_value_skat, weight,
  q_value`, sorted by `p_value` ascending.
  - There is no `effect_size` column: SKAT-O is a variance-component test
    and does not produce one.
  - `q_value` is a real Benjamini-Hochberg FDR correction (`bh_fdr`) over
    `p_value`, not the old mislabeled Bonferroni value.

Not yet ported from the original pipeline (future work):
- kinship/familial null models (the current null model assumes unrelated
  samples)
- sliding-window region grouping
- the Fisher/Barnard contingency-table testing arm (`method="fisher"` is
  declared in `StatsConfig` but not implemented)

## Backend Boundary
The contracts above are backend-neutral.
Scheduler-specific concerns (Slurm, AWS Batch, etc.) must not change stage I/O schemas.
