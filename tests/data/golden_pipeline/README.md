# Golden Pipeline Fixture (PH2-013)

A frozen, versioned end-to-end fixture for the `bed_to_training -> train ->
predict -> aggregate -> stats` pipeline, used by the golden regression tests
(`tests/test_golden_fixture_integrity.py`, `tests/test_golden_predict.py`,
`tests/test_golden_aggregate_stats.py`) and by `tests/conftest.py`'s
`golden_*` fixtures.

## Provenance

`genome.fa` and `training/` were fetched once from Ensembl REST via:

```
python scripts/fetch_golden_genome_fixture.py \
    --chrom chr22 --start 20000000 --end 20100000 \
    --out tests/data/golden_pipeline
```

- **FASTA**: GRCh38 `chr22:20,000,000-20,100,000` (100,000 bp, no `N`
  bases), fetched via `GET /sequence/region/human/22:...` (Ensembl REST).
- **Training regions**: Ensembl Regulatory Build features
  (`enhancer` / `CTCF_binding_site` / `promoter`) overlapping the same
  region, fetched via `GET /overlap/region/human/22:...?feature=regulatory`.
  32 regions, 3 features — small enough that no trimming was needed to
  bound 2-epoch CPU training time.
- `genome.fa.fai` was built via `pyfaidx.Faidx` and is committed alongside
  `genome.fa` so pyfaidx never writes into the working tree during test
  runs (which would otherwise dirty `git status` on every run).

This replaced a prior design where `tests/conftest.py`'s `real_fasta` /
`real_training_dir` fixtures fetched this same data live from Ensembl REST
on every test session — a genuine flake risk (a live 503 was observed from
that code) and a hard requirement for network access to run the suite at
all.

## Coordinates are LOCAL to this fixture

**Position 1 in `variants.tsv`, or position 0 in `training/training_regions.bed`,
is NOT a real chr22 coordinate.** All BED/variant coordinates in this
fixture are relative to the start of `genome.fa` — position 0 (BED) / VCF
POS 1 (1-based) corresponds to real GRCh38 `chr22:20,000,000`. This matches
the convention the old live-fetch `conftest.py` fixtures already used; it
is preserved here so nothing downstream needs to know the difference.

## Inputs (hand-maintained) vs. outputs (regenerated)

**Hand-maintained inputs** — changed only by deliberate edit, never by
re-running a generator blindly:
- `genome.fa`, `genome.fa.fai`, `training/training_regions.bed`,
  `training/features.txt` — via `scripts/fetch_golden_genome_fixture.py`
  (rerun only if the frozen region needs to change).
- `variants.tsv` — via `scripts/build_golden_variants.py`. 18 variants (14
  SNV, 2 insertion, 2 deletion) at training-region midpoints, each with a
  500bp flank clear of both fixture edges. Reference alleles are read
  directly from `genome.fa` via pyfaidx (never hand-typed), so the
  pipeline's own reference-allele-mismatch sanity-check warning never
  fires against this fixture (see `tests/test_golden_predict.py`, which
  asserts exactly that via `caplog`).
- `predictions.tsv` — via `scripts/build_golden_predictions.py`. Values
  are an explicit deterministic formula, **not** run through an actual
  model (an untrained/random model produces near-identical values across
  variants/features, which would hide column-mixup bugs). Column layout
  matches `VariantEffectPredictor.predict_variants`'s real emission order:
  the variant columns followed by one `ref_pred_<i>`/`alt_pred_<i>`/
  `delta_<i>` triple per model feature (not grouped by kind).
- `sample_ids.txt`, `genotypes.tsv`, `phenotype.tsv` — hand-authored
  (deterministic, not from a generator script). **Synthetic**: these are
  not real cohort genotypes/phenotypes. Real cohort statistics are
  PH2-020 Phase 2, deliberately out of scope here (see
  `docs/2026-08-13-session-handover.md` and `PROJECT_PLAN.md`).

**Regenerated outputs** — never hand-edited, only written by
`scripts/regenerate_golden_fixtures.py`:
- `expected/aggregate_scores.tsv`, `expected/aggregate_genotypes.tsv`,
  `expected/stats_results.tsv` — the `aggregate` -> `stats` stages' real
  output, run in-process via `LocalRunner` with the R-backed `RpyBackend`
  swapped for `tests/golden_utils.DeterministicSkatBackend` (a pure,
  input-sensitive stand-in — not a statistical test — so this fixture
  needs no R/rpy2 installed).

To regenerate (after a deliberate, understood pipeline change):

```
python scripts/regenerate_golden_fixtures.py            # writes expected/*.tsv
python scripts/regenerate_golden_fixtures.py --check     # diffs instead of writing
```

**Before committing regenerated `expected/*.tsv` files, hand-verify the
arithmetic on at least 2 rows** against `predictions.tsv` — a wrong golden
file would silently bless a bug. This is not optional; see the git history
for a worked example.

## Two pinned-not-fixed bugs

This fixture deliberately locks in two known-real, not-yet-fixed pipeline
behaviors rather than hiding them, so a future change to either is a
visible, deliberate diff instead of a silent regression:

1. **Missing genotype call recodes to dosage 9** (not NaN/excluded):
   `genotypes.tsv` has exactly one `./.` call, at variant
   `chr22:13492:G:T` (GENE_A, variant index 5 in `variants.tsv`), sample
   `S05`.
2. **Phenotype-absent samples are silently dropped**: `phenotype.tsv`
   omits sample `S09`, which IS present in `genotypes.tsv`/
   `sample_ids.txt`. The pipeline does not warn or error — it silently
   analyzes on the 8 remaining samples.

See `PROJECT_PLAN.md`'s PH2-013 entry for follow-up tracking.
