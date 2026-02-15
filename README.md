## Current state
This code has not been tested since the refactoring. This is going to happen, when a large chunk of the code has been moved.

# PHDeep - Pulmonary Hypertension Deep Learning pipeline
This project is a refactored version of my PhD project, based on the thesis "Deep Learning of regulatory sequence variation in Pulmonary Arterial Hypertension". It is still strongly under development.

# Project structure
This project is modularised to be as adjustable as possible, while still (hopefully) easy to understand. Most module names should be self-explanatory. Something to point out are prediction and post_prediction. The former holding functions to actually run across either the reference genome or through vcf files, fetching the reference around the variants on the go. For this, the data_loading module is used to load the files. 
Post_prediction is for the handling of predicted epigenetic marks. Here, you will find the functions that compute the scores and make them accessible to the statistical/association testing. 

# Future work
These modules will be connected in Nextflow Pipelines, so that it is possible to run the pipeline from training-set generation over model training to prediction and post_prediction. 

## Development setup with uv

1. Install `uv`:

```bash
python3 -m pip install --user uv
```

2. Create and sync the environment:

```bash
uv venv .venv
uv sync --extra dev
```

3. Activate:

```bash
source .venv/bin/activate
```

## ENCODE fixture for bed_to_training

Generate a tiny reproducible fixture from ENCODE:

```bash
.venv/bin/python scripts/fetch_encode_dnase_fixture.py
```

Run the converter against the fixture:

```bash
.venv/bin/python src/data_processing/bed_to_training.py \
  -m tests/data/encode_dnase_fixture/metadata.tsv \
  -i tests/data/encode_dnase_fixture/input_bed_files \
  -o tests/data/encode_dnase_fixture/output/training_regions.bed \
  -e tests/data/encode_dnase_fixture/output/read_errors.txt \
  -f tests/data/encode_dnase_fixture/output/features.txt \
  -a GRCh38
```

Run tests:

```bash
.venv/bin/python -m pytest -q
```

Validate pipeline config contract:

```bash
.venv/bin/python main.py --validate-config config/pipeline.example.json
```
