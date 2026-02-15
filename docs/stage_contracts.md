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

Outputs:
- table with minimum columns:
  - `chromosome`
  - `start`
  - `end`
  - `prediction`

## Stage: `aggregate`
Description:
Aggregates per-position predictions to per-feature/per-region scores.

Inputs:
- one or more prediction tables from `predict`

Outputs:
- score table with minimum columns:
  - `chromosome`
  - `start`
  - `end`
  - `score`

## Stage: `stats`
Description:
Consumes aggregated score tables for association testing.

Inputs:
- aggregated score table
- method config (`skat_o`, `fisher`, or `none`)

Outputs:
- result table with minimum columns:
  - `feature_id`
  - `p_value`
  - `q_value`
  - `effect_size`

## Backend Boundary
The contracts above are backend-neutral.
Scheduler-specific concerns (Slurm, AWS Batch, etc.) must not change stage I/O schemas.
