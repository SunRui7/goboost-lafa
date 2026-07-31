# GOBoost FunctionBench / LAFA submission

Public source repository:
[SunRui7/goboost-lafa](https://github.com/SunRui7/goboost-lafa).

## Artifact to submit

FunctionBench evaluates a runnable container image rather than a one-time prediction file.
Build this repository's `Dockerfile`, run the smoke test below, and publish the resulting
image to a public Docker Hub repository. Give the LAFA maintainers both the public image tag
and the public source repository URL.

## Interface

Required arguments:

- `--query_file`, `-q`: query sequences in FASTA or FASTA.gz format
- `--graph`: current Gene Ontology graph in OBO or OBO.gz format
- `--output_baseline`, `--output_file`, `-o`: output TSV or TSV.gz path

Compatibility arguments accepted by the current LAFA backend:

- `--annot_file`, `-a`
- `--train_sequences`
- `--train_taxonomy`
- `--num_threads`

Optional GOBoost arguments:

- `--structure_dir`: directory containing extracted AlphaFold/PDB files
- `--esm_model`: mounted ESM-1b checkpoint (default
  `/app/weights/esm1b_t33_650M_UR50S.pt`)
- `--device`: `auto` (default), `cpu`, `cuda`, or a CUDA device
- `--contact_threshold`: ESM contact edge threshold, default `0.5`
- `--min_score`: minimum emitted score, default `0.01`
- `--n_terms`: maximum terms per protein, default `1500`; `0` disables the cap

Output is a headerless, tab-separated file with exactly three columns:

```text
P12345\tGO:0005737\t0.123456
P12345\tGO:0016020\t0.098765
```

Query identifiers follow the LAFA baseline convention: for UniProt headers such as
`sp|P12345|NAME`, the emitted identifier is `P12345`.

## Reproducibility and data

- The nine published GOBoost `All`, `Head`, and `Tail` checkpoints are downloaded from
  [Zenodo record 14048928](https://zenodo.org/records/14048928) during the image build.
- The Docker build verifies `Model.zip` with the Zenodo-published MD5 checksum
  `bbdf045c13db349acab8b7777230184d`.
- ESM-1b is 7.29 GiB and is therefore an external read-only mount, following LAFA's
  guidance for data larger than 5 GB. Download both files from the official fair-esm host:
  `https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt` and
  `https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt`.
  The local model MD5 is
  `ba8914bc3358cae2254ebc8874ee67f6`; the regression MD5 is
  `e7fe626dfd516fb6824bd1d30192bdb1`.
- With the ESM volume mounted, runtime inference requires no network access.
- The container can run from FASTA alone. AlphaFoldDB v6 structures are optional and can be
  mounted separately to avoid enlarging the image.
- Sequences are limited to 1,000 residues, matching the published GOBoost predictor.
- GOBoost is a frozen model: timepoint training files are not used to update its parameters.

## Build and smoke test

```bash
docker build -t aubrey7/goboost-lafa:v1 .

mkdir -p output
docker run --rm --gpus all \
  -v /absolute/path/to/lafa-release:/app/data:ro \
  -v /absolute/path/to/weights:/app/weights:ro \
  -v "$PWD/output:/app/output:rw" \
  aubrey7/goboost-lafa:v1 \
  --query_file /app/data/test_sequences.fasta \
  --train_sequences /app/data/train_sequences.fasta \
  --annot_file /app/data/train_terms.tsv \
  --graph /app/data/go-basic.obo \
  --output_baseline /app/output/goboost_predictions.tsv.gz \
  --num_threads 16

python scripts/validate_lafa_output.py output/goboost_predictions.tsv.gz
```

For an optional structure-backed run, mount the directory and add `--structure_dir`:

```bash
-v /absolute/path/to/AlphaFoldDB_v6:/app/structures:ro
--structure_dir /app/structures
```

## Publish

```bash
docker login
docker push aubrey7/goboost-lafa:v1
docker pull aubrey7/goboost-lafa:v1
```

Make the image public, add the example run command to its Docker Hub description, and link
the source repository. Use [DOCKERHUB_OVERVIEW.md](DOCKERHUB_OVERVIEW.md) as the Docker Hub
description template. The FunctionBench **Submit Your Method** link currently points to the
LAFA container guide; contact `idoerg@iastate.edu` with the public image tag and public
source URL once the anonymous pull test succeeds.
