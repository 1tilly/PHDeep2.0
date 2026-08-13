# PHDeep2.0 Project Plan

Last updated: 2026-02-15

## Goal
Reach functional parity with the original `PhDeep` feature set while making `PHDeep2.0` infrastructure-independent:
- No hardcoded Cambridge HPC paths.
- No Slurm-specific logic inside core modules.
- Support local execution first, with optional adapters for Slurm/AWS/other schedulers.

## Product Principles
1. Core library is scheduler-agnostic.
2. All paths are config/env-driven.
3. One primary model stack in core (recommend PyTorch).
4. Reproducibility and tests are mandatory for every merged feature.
5. Workflow orchestration is an outer layer, not embedded in business logic.

## Progress Log
### 2026-02-15
- Added packaging baseline with `pyproject.toml`.
- Added baseline runtime dependencies to `requirements.txt`.
- Added package init files under `src/` and submodules for import/package discovery.
- Replaced empty `main.py` with a basic CLI entrypoint and path resolution smoke command.
- Replaced empty `config/paths.py` with env-driven path config (`PHDEEP_*` variables).
- Fixed major loader/runtime blockers:
  - `src/data_loading/bed_loader.py`: constructor state, logging API bug, argument mismatch, comparison bugs.
  - `src/data_processing/bed_to_training.py`: unresolved function calls and nonfunctional CLI flow.
  - `src/data_processing/vcf_processing.py`: invalid inheritance and string append bug.
  - `src/data_loading/vcf_loader.py`: removed hardcoded `R_HOME` and bcftools binary paths; added optional `rpy2` handling.
  - `src/models/base_model/architecture.py`: abstract/classmethod and output-shape computation fix.
  - `src/models/blainville/architecture.py` and `src/models/jellyfish/architecture.py`: fixed base import path and jellyfish kernel variable bugs.
- Verified syntax with `python -m compileall -q .` (pass).
- Added explicit stage contract document: `docs/stage_contracts.md`.
- Added typed pipeline config schema and validator: `config/pipeline_config.py`.
- Added executable example config: `config/pipeline.example.json`.
- Added tests for contract/schema validation: `tests/test_pipeline_config.py`.
- Added backend-neutral runner interface with `LocalRunner`: `src/workflow/runners.py`.
- Wired `main.py --run-config` to execute stage graph from typed config.
- Added runner tests: `tests/test_workflow_runner.py`.
- Hardened `bed_to_training` by replacing deprecated dataframe append flow with concat.
- Added portability guardrail test for hardcoded site paths: `tests/test_portability_guardrails.py`.

### 2026-08-13
- Rewrote `src/post_prediction/aggregation.py` (PH2-019): added
  `build_variant_weights_table` (per-variant weights, not collapsed to
  one row per group) and `build_genotype_matrix`; fixed the doubled-label
  bug in `aggregate_variant_scores`; removed the dead/wrong-shape
  `build_skat_input`.
- Rewrote `src/statistical_testing/skat_o_test.py` (PH2-020): `run_skat_o`
  and `run_skat_o_from_feather` now take real genotype/phenotype inputs
  and a pluggable `SkatBackend` (module import never requires rpy2); added
  `bh_fdr` for real Benjamini-Hochberg q-values.
- Rewired `config/pipeline_config.py` and `src/workflow/runners.py` to the
  corrected `aggregate`/`stats` contracts (new `AggregateConfig`/
  `StatsConfig` fields, `_run_aggregate`/`_run_stats` updated to match,
  `group_col`/`weight_col`/`min_variants` now actually threaded through to
  `run_skat_o_from_feather`), and corrected the `predict`/`aggregate`/
  `stats` sections of `docs/stage_contracts.md`.
- Phase 2 remains deferred future work: real genotype extraction from a
  cohort's actual VCF/BCF data (today's `build_genotype_matrix` only
  recodes genotype columns already present on the input DataFrame), cohort
  phenotype/covariate ingestion, and familial/kinship null models.

## Milestones

### M0 - Baseline and Guardrails (Week 1)
Scope:
- Stabilize imports and runtime blockers.
- Add packaging and dependency management.
- Add lint/type/test tooling skeleton.

Exit criteria:
- `python -m compileall` passes.
- Import smoke tests pass for non-optional modules.
- Local setup from clean environment is documented and works.

### M1 - Portability Foundation (Week 1-2)
Scope:
- Introduce typed configuration (`paths`, `runtime`, `execution backend`).
- Remove hardcoded absolute paths from code.
- Define backend interface: `local`, `slurm`, `aws_batch` adapters.

Exit criteria:
- No committed code contains hardcoded `/home/...`, `/rds/...`, or site-specific binaries.
- Same command can run on local backend and mocked batch backend via config only.

### M2 - Data Pipeline Parity (Week 2-4)
Scope:
- BED metadata -> training BED/features pipeline parity.
- GFF/GeneHancer parsing hardening.
- Genome and variant sequence generation hardening.

Exit criteria:
- Golden dataset outputs match agreed reference outputs.
- Data processing steps have unit tests and integration tests.

### M3 - Model Training Parity (Week 4-6)
Scope:
- Implement one canonical DeepSEA-family model end-to-end first.
- Training loop, evaluation metrics, checkpointing, deterministic seeds.
- Experiment config schema and CLI.

Exit criteria:
- Train/eval works locally on fixture dataset.
- Artifacts and metrics are reproducible with fixed seed.

### M4 - Prediction and Post-Prediction Parity (Week 6-8)
Scope:
- Reference and variant prediction pipelines.
- Aggregation and scoring outputs for association testing handoff.
- Optional BCF/R integration moved behind adapters.

Exit criteria:
- End-to-end mini run: predict -> aggregate -> stats input export.

### M5 - Workflow Adapters and Deployment Paths (Week 8-9)
Scope:
- Backend runners for local, Slurm, and AWS Batch.
- Parameterized templates, no project-specific paths.

Exit criteria:
- Same pipeline config runs unchanged on at least local + one batch backend.

### M6 - Hardening, CI, and Release (Week 9-11)
Scope:
- Regression suite against golden fixtures.
- CI gates: lint, unit/integration tests, packaging checks.
- Migration docs from original project.

Exit criteria:
- CI is green.
- Release candidate tagged with upgrade notes.

## Issue Backlog

Status legend: `todo`, `in_progress`, `blocked`, `done`

| ID | Milestone | Title | Priority | Depends on | Status | Acceptance criteria |
|---|---|---|---|---|---|---|
| PH2-001 | M0 | Create package skeleton and `pyproject.toml` | P0 | - | done | Editable install works; basic CLI entrypoint exists |
| PH2-002 | M0 | Add pinned dependencies and optional extras (`slurm`, `aws`, `r`) | P0 | PH2-001 | in_progress | Clean install passes import smoke tests |
| PH2-003 | M0 | Fix broken imports/signatures in loaders/models | P0 | PH2-001 | in_progress | Current runtime blockers resolved |
| PH2-004 | M0 | Add pre-commit hooks (ruff/black/mypy/pytest) | P1 | PH2-001 | todo | Hooks run locally and in CI |
| PH2-005 | M1 | Implement typed settings module (`paths`, `env`, `backend`) | P0 | PH2-001 | done | All path/runtime values loaded from config/env |
| PH2-006 | M1 | Remove hardcoded absolute paths from source tree | P0 | PH2-005 | done | `rg '/home/|/rds/'` returns no code hits |
| PH2-007 | M1 | Define execution backend interface | P0 | PH2-005 | done | `LocalRunner` implemented with contract tests |
| PH2-008 | M1 | Add backend adapters (`SlurmRunner`, `AwsBatchRunner`) as optional modules | P1 | PH2-007 | todo | Adapters can be imported without affecting core |
| PH2-009 | M2 | Repair and finalize BED parser API | P0 | PH2-003, PH2-005 | in_progress | BED parser passes unit tests on fixtures |
| PH2-010 | M2 | Rebuild `bed_to_training` workflow using parser API | P0 | PH2-009 | in_progress | Feature list + training bed outputs reproducible |
| PH2-011 | M2 | Harden genome/variant sequence generation (`vcf_processing`) | P0 | PH2-003 | todo | Sequence mutation cases (SNV/INS/DEL) validated |
| PH2-012 | M2 | Harden GeneHancer/GFF parser behavior and errors | P1 | PH2-003 | todo | Unsupported paths raise explicit typed errors |
| PH2-013 | M2 | Build golden fixture dataset and expected outputs | P0 | PH2-009, PH2-010, PH2-011 | todo | Fixtures versioned; regression tests pass |
| PH2-014 | M3 | Fix model module import structure and package paths | P0 | PH2-003, PH2-001 | todo | Model modules import by package path only |
| PH2-015 | M3 | Correct channel-size calculation and model shape tests | P0 | PH2-014 | todo | Forward pass succeeds with shape assertions |
| PH2-016 | M3 | Implement canonical training loop CLI and checkpointing | P0 | PH2-014, PH2-015, PH2-005 | todo | `train` command runs end-to-end locally |
| PH2-017 | M3 | Add deterministic seed and reproducibility checks | P1 | PH2-016 | todo | Repeat runs within expected variance |
| PH2-018 | M4 | Implement prediction pipeline CLI (`reference`, `variant`) | P0 | PH2-011, PH2-016 | todo | Prediction artifacts produced for fixtures |
| PH2-019 | M4 | Implement post-prediction aggregation and scoring | P0 | PH2-018 | done | Aggregated outputs match schema and tests |
| PH2-020 | M4 | Define stats handoff contract and export format | P1 | PH2-019 | done | Contract doc + integration tests complete |
| PH2-021 | M4 | Add optional R/BCF adapter wrappers | P2 | PH2-020, PH2-007 | todo | Optional adapters enabled by extras only |
| PH2-022 | M5 | Implement local workflow runner over stage graph | P0 | PH2-007, PH2-010, PH2-016, PH2-019 | in_progress | Full mini pipeline runnable from one command |
| PH2-023 | M5 | Implement Slurm workflow adapter using same stage graph | P1 | PH2-022, PH2-008 | todo | Same graph executes via Slurm config |
| PH2-024 | M5 | Implement AWS Batch workflow adapter using same stage graph | P1 | PH2-022, PH2-008 | todo | Same graph executes via AWS Batch config |
| PH2-025 | M6 | CI pipeline with lint/type/tests/build matrix | P0 | PH2-004, PH2-013, PH2-016, PH2-019 | todo | CI required checks enforced on main |
| PH2-026 | M6 | Migration guide from original `PhDeep` | P1 | PH2-010, PH2-016, PH2-019 | todo | Docs include command mapping + caveats |
| PH2-027 | M6 | Release candidate checklist and tagging process | P1 | PH2-025, PH2-026 | todo | v0 parity RC tagged with release notes |

## Dependency Graph

```text
PH2-001 -> PH2-002 -> PH2-003
PH2-001 -> PH2-004

PH2-001 -> PH2-005 -> PH2-006
PH2-005 -> PH2-007 -> PH2-008

PH2-003 + PH2-005 -> PH2-009 -> PH2-010
PH2-003 -> PH2-011
PH2-003 -> PH2-012
PH2-009 + PH2-010 + PH2-011 -> PH2-013

PH2-001 + PH2-003 -> PH2-014 -> PH2-015 -> PH2-016 -> PH2-017

PH2-011 + PH2-016 -> PH2-018 -> PH2-019 -> PH2-020 -> PH2-021

PH2-007 + PH2-010 + PH2-016 + PH2-019 -> PH2-022
PH2-022 + PH2-008 -> PH2-023
PH2-022 + PH2-008 -> PH2-024

PH2-004 + PH2-013 + PH2-016 + PH2-019 -> PH2-025
PH2-010 + PH2-016 + PH2-019 -> PH2-026
PH2-025 + PH2-026 -> PH2-027
```

## Risks and Mitigations
1. Risk: Scope creep from trying to port all legacy scripts.
Mitigation: Only parity against agreed feature matrix; keep `scratch` out of critical path.

2. Risk: Reintroducing HPC coupling via convenience shortcuts.
Mitigation: CI check for forbidden hardcoded path patterns and backend imports in core modules.

3. Risk: Mixed frameworks increase complexity and maintenance cost.
Mitigation: Standardize core training/prediction on one framework; keep others behind adapters only if needed.

4. Risk: Regression uncertainty due to limited legacy tests.
Mitigation: Golden fixtures and contract tests for each pipeline stage.

## Definition of Done (Parity Release)
1. Core stages run locally from clean environment with documented commands.
2. No infrastructure-specific paths in core code.
3. Workflow can target local plus at least one batch scheduler by config.
4. Golden fixture outputs and integration tests pass in CI.
5. Migration doc covers equivalent commands and known differences from original `PhDeep`.
