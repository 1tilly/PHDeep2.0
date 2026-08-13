# PHDeep2.0 — Deep Learning of Regulatory Sequence Variation in Pulmonary Arterial Hypertension

A portable, modular PyTorch pipeline for predicting the regulatory impact of
non-coding genetic variants, developed as part of a PhD thesis on Pulmonary
Arterial Hypertension (PAH).

The pipeline covers the full workflow from ENCODE metadata to variant effect
scores ready for SKAT-O association testing:

```
ENCODE metadata + BED files
        ↓  bed_to_training
Training regions (BED) + genome FASTA
        ↓  GenomicsDataset + Trainer
Trained model (checkpoint)
        ↓  VariantEffectPredictor
Per-variant Δ scores (alt − ref)
        ↓  aggregation
Per-variant weights (grouped by gene/region) → SKAT-O
```

## Model zoo

| Name | Channels | Conv blocks | Extra features |
|---|---|---|---|
| `DeepSEA` | 320 → 480 → 960 | 3 × 1 conv | Original architecture (Zhou & Troyanskaya, 2015) |
| `DeeperDeepSEA` | 320 → 480 → 960 | 3 × 2 conv | Doubled depth, BatchNorm after pooling |
| `BlainvilleDeepSEA` | 320 → 480 → 960 | 3 + 3 + 2 conv | Tripled first block |
| `BottleNoseDeepSEA` | 160 → 320 → 720 → 1080 | 4 blocks | Wider channels, bottleneck expansion |
| `JellyFishDeepSEA` | configurable | configurable | Parameterised depth/width, WandB-sweep friendly |

All models share the `AbstractCNN` base class and accept `(sequence_length, n_targets)`.

## Installation

Requires Python ≥ 3.10 and `bedtools`/`tabix` on PATH.

```bash
# 1. Create environment
python -m venv .venv && source .venv/bin/activate

# 2. Install with model and dev extras
pip install -e ".[models,dev]"

# Optional: R integration for SKAT-O
pip install -e ".[stats]"
```

With `uv` (faster):

```bash
uv venv .venv
uv sync --extra dev --extra models
source .venv/bin/activate
```

## Quick start

### 1. Build the training BED from ENCODE metadata

```bash
python src/data_processing/bed_to_training.py \
  -m data/encode_metadata.tsv \
  -i data/bed_files/ \
  -o output/training_regions.bed \
  -e output/read_errors.txt \
  -f output/features.txt \
  -a GRCh38
```

A reproducible fixture can be fetched from ENCODE for testing:

```bash
python scripts/fetch_encode_dnase_fixture.py
```

### 2. Train a model

```python
from torch.utils.data import DataLoader, random_split

from src.data_loading.genomics_dataset import GenomicsDataset
from src.models.deepsea.architecture import DeepSEA
from src.training.trainer import Trainer, get_criterion, get_optimizer, set_seed

set_seed(42)

with open("output/features.txt") as f:
    feature_list = [l.strip() for l in f if l.strip()]

dataset = GenomicsDataset(
    bed_path="output/training_regions.bed",
    genome_fasta="data/GRCh38.fa",
    feature_list=feature_list,
    seq_len=1000,
)
n_val = int(0.1 * len(dataset))
train_ds, val_ds = random_split(dataset, [len(dataset) - n_val, n_val])

train_loader = DataLoader(train_ds, batch_size=128, shuffle=True, num_workers=4)
val_loader   = DataLoader(val_ds,   batch_size=256, shuffle=False, num_workers=4)

model   = DeepSEA(sequence_length=1000, n_targets=len(feature_list))
opt     = get_optimizer(model, lr=0.01)
loss_fn = get_criterion()

trainer = Trainer(model, opt, loss_fn, device="cuda", checkpoint_dir="checkpoints/")
history = trainer.fit(train_loader, val_loader, n_epochs=100)
```

### 3. Predict variant effects

```python
import torch
from src.models.deepsea.architecture import DeepSEA
from src.prediction.predict import VariantEffectPredictor

model = DeepSEA(sequence_length=1000, n_targets=len(feature_list))
predictor = VariantEffectPredictor(
    model,
    genome_fasta="data/GRCh38.fa",
    seq_len=1000,
    checkpoint_path="checkpoints/best_model.pt",
)
delta_df = predictor.predict_variants(var_df)
```

### 4. Validate pipeline config and run stages

```bash
python main.py --validate-config config/pipeline.example.json
python main.py --run-config config/pipeline.example.json
```

### 5. Run tests

```bash
pytest -q
```

## Project structure

```
src/
  data_loading/
    bed_loader.py          ENCODE BED metadata parser
    genome_loader.py       FASTA sequence access and one-hot encoding
    genomics_dataset.py    PyTorch Dataset for multi-label training
    vcf_loader.py          BCFtools VCF/BCF query wrapper
    gff_loader.py          GeneHancer GFF parser
    ensembl_api.py         Ensembl REST API client
  data_processing/
    bed_to_training.py     Build labelled training BED from ENCODE metadata
    vcf_processing.py      Variant sequence extraction (SNV/INS/DEL)
  models/
    base_model/            AbstractCNN base class
    deepsea/               Canonical DeepSEA (Zhou & Troyanskaya 2015)
    deeper_deepsea/        Doubled-depth DeepSEA
    blainville/            Tripled first-block variant
    bottlenose/            Wide-channel variant (1080 channels)
    jellyfish/             Configurable depth/width model
    bert_models/           BERT-based sequence model (experimental)
  training/
    trainer.py             Training loop, early stopping, checkpointing
  prediction/
    predict.py             Reference and variant effect prediction
  post_prediction/
    aggregation.py         Per-variant weight table + genotype matrix aggregation
  statistical_testing/
    skat_o_test.py         SKAT-O wrapper (requires R extra)
  workflow/
    runners.py             Backend-neutral pipeline runner (local/Slurm/AWS)
config/
  pipeline_config.py       Typed pipeline config schema
  pipeline.example.json    Example pipeline configuration
tests/                     Pytest test suite
scripts/                   Utility scripts (fixture fetching, etc.)
```

## Background

This project is a re-implementation of the computational pipeline developed
for the thesis *"Deep Learning of regulatory sequence variation in Pulmonary
Arterial Hypertension"*. The original pipeline ran on the Cambridge HPC
cluster with hardcoded paths and Slurm dependencies. PHDeep2.0 removes all
infrastructure coupling: all paths are config- or environment-driven, and
the scheduler is an optional adapter layer.

The models follow the DeepSEA family (Zhou & Troyanskaya, Nature Methods
2015) — convolutional networks trained to predict chromatin accessibility,
transcription factor binding, and histone marks from raw DNA sequence.
Variant effect scores are computed as the difference in predicted activity
between reference and alternate alleles, then turned into a per-variant
weights table (grouped by whatever `group_col` the caller supplies — e.g.
gene symbol) for burden/variance-component testing via SKAT-O.

## References

- Zhou, J. & Troyanskaya, O.G. (2015). Predicting effects of noncoding
  variants with deep learning-based sequence model. *Nature Methods*, 12,
  931–934. https://doi.org/10.1038/nmeth.3547
- Kelley, D.R. et al. (2018). Sequential regulatory activity prediction
  across chromosomes with convolutional neural networks. *Genome Research*.
  https://doi.org/10.1101/gr.227819.117

## License

MIT — see [LICENSE](LICENSE).
