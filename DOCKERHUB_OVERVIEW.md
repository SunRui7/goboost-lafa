# GOBoost-LAFA

LAFA-compatible container for GOBoost protein function prediction.

- Source: https://github.com/SunRui7/goboost-lafa
- Original method: https://github.com/Cao-Labs/GOBoost
- Container interface: https://github.com/anphan0828/LAFA_container_guide

## Required external weights

ESM-1b is larger than 5 GB and is intentionally supplied as an external,
read-only volume. Download both files and keep their names unchanged:

- https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
- https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt

Mount the directory containing both files at `/app/weights`.

## Run

```bash
docker run --rm --gpus all \
  -v /path/to/lafa-release:/app/data:ro \
  -v /path/to/esm-weights:/app/weights:ro \
  -v /path/to/output:/app/output:rw \
  aubrey7/goboost-lafa:v1 \
  --query_file /app/data/test_sequences.fasta \
  --train_sequences /app/data/train_sequences.fasta \
  --annot_file /app/data/train_terms.tsv \
  --graph /app/data/go-basic.obo \
  --output_file /app/output/goboost_predictions.tsv.gz \
  --num_threads 16
```

The output is a gzip-compressed, headerless, three-column TSV containing
`Query_ID`, `GO_Term`, and `Score`.

## Optional AlphaFoldDB structures

The container can run from FASTA alone using ESM-1b embeddings and predicted
residue contacts. To use AlphaFoldDB v6 structures, mount the extracted PDB
directory read-only at `/app/structures` and add:

```text
--structure_dir /app/structures
```

Supported structure names include
`AF-<EntryID>-F1-model_v6.pdb.gz`, `AF-<EntryID>-F1-model_v4.pdb.gz`,
`<EntryID>.pdb.gz`, and `<EntryID>.pdb`.

## Output contract

- Headerless TSV or TSV.gz
- Exactly three columns: query ID, GO term, score
- Scores are in `[0, 1]`
- UniProt FASTA headers such as `sp|P12345|NAME` emit the ID `P12345`
