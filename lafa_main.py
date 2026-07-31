#!/usr/bin/env python3
"""FunctionBench/LAFA entry point for the published GOBoost checkpoints."""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import numpy as np


ASPECTS = ("mf", "bp", "cc")
OUTPUT_DIMS = {"bp": (1943, 300, 1643), "mf": (489, 100, 389), "cc": (320, 50, 270)}
ROOT_TERMS = {"bp": "GO:0008150", "mf": "GO:0003674", "cc": "GO:0005575"}
AA_ALPHABET = "ARNDCQEGHILKMFPSTWYVX"
AA_TO_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GOBoost on a FASTA file and emit headerless LAFA/CAFA TSV predictions."
    )
    parser.add_argument("--query_file", "-q", required=True, help="Query protein FASTA (optionally gzipped)")
    parser.add_argument("--graph", required=True, help="GO ontology in OBO format")
    parser.add_argument(
        "--output_baseline", "--output_file", "-o", dest="output_file", required=True,
        help="Headerless three-column prediction TSV (optionally .gz)",
    )
    # Accepted for compatibility with the standard LAFA invocation. GOBoost is
    # a frozen model and does not retrain at each timepoint.
    parser.add_argument("--annot_file", "-a", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--train_sequences", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--train_taxonomy", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--model_dir", default="/app/Model", help="Directory containing GOBoost .pt files")
    parser.add_argument(
        "--esm_model", default=os.environ.get("GOBOOST_ESM_MODEL", "/app/weights/esm1b_t33_650M_UR50S.pt"),
        help="Mounted local ESM-1b checkpoint; the contact-regression file must be beside it",
    )
    parser.add_argument(
        "--structure_dir", default=None,
        help="Optional AlphaFold/PDB directory; ESM contact predictions are used when a structure is absent",
    )
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, or a CUDA device such as cuda:0")
    parser.add_argument("--contact_threshold", type=float, default=0.5)
    parser.add_argument("--min_score", type=float, default=0.01)
    parser.add_argument("--n_terms", type=int, default=1500, help="Maximum predictions per query; 0 means unlimited")
    parser.add_argument("--num_threads", type=int, default=8)
    parser.add_argument("--max_length", type=int, default=1000, choices=[1000])
    args = parser.parse_args(argv)
    if not 0 <= args.contact_threshold <= 1:
        parser.error("--contact_threshold must be in [0, 1]")
    if not 0 <= args.min_score <= 1:
        parser.error("--min_score must be in [0, 1]")
    if args.n_terms < 0 or args.num_threads < 1:
        parser.error("--n_terms must be non-negative and --num_threads must be positive")
    return args


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def entry_id(header: str) -> str:
    token = header.split(None, 1)[0]
    fields = token.split("|")
    return fields[1] if len(fields) > 1 and fields[1] else token


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    sequence: List[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield entry_id(header), "".join(sequence).upper()
                header = line[1:].strip()
                sequence = []
            elif header is None:
                raise ValueError(f"{path}:{line_number}: sequence found before the first FASTA header")
            else:
                sequence.append("".join(line.split()))
    if header is not None:
        yield entry_id(header), "".join(sequence).upper()


def load_term_vocabulary(path: Path) -> Dict[str, List[str]]:
    labels: Dict[str, List[str]] = {}
    with path.open(encoding="utf-8") as handle:
        lines = [next(handle, "").rstrip("\n") for _ in range(10)]
    for aspect, index in (("mf", 1), ("bp", 5), ("cc", 9)):
        labels[aspect] = [term for term in lines[index].split("\t") if term]
        expected = OUTPUT_DIMS[aspect][0]
        if len(labels[aspect]) != expected:
            raise ValueError(f"Expected {expected} {aspect.upper()} terms, found {len(labels[aspect])}")
    return labels


def load_ontology(path: Path) -> Tuple[set[str], Dict[str, set[str]]]:
    active: set[str] = set()
    parents: Dict[str, set[str]] = defaultdict(set)
    term_id = None
    obsolete = False
    term_parents: set[str] = set()

    def commit() -> None:
        if term_id and not obsolete:
            active.add(term_id)
            parents[term_id].update(term_parents)

    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line == "[Term]":
                commit()
                term_id, obsolete, term_parents = None, False, set()
            elif line.startswith("["):
                commit()
                term_id, obsolete, term_parents = None, False, set()
            elif line.startswith("id: GO:"):
                term_id = line.split()[1]
            elif line == "is_obsolete: true":
                obsolete = True
            elif line.startswith("is_a: GO:"):
                term_parents.add(line.split()[1])
            elif line.startswith("relationship: part_of GO:"):
                term_parents.add(line.split()[2])
    commit()
    if not active:
        raise ValueError(f"No active GO terms found in {path}")
    return active, parents


def ancestors(term: str, parents: Mapping[str, set[str]], memo: Dict[str, set[str]]) -> set[str]:
    if term in memo:
        return memo[term]
    result: set[str] = set()
    stack = list(parents.get(term, ()))
    while stack:
        parent = stack.pop()
        if parent not in result:
            result.add(parent)
            stack.extend(parents.get(parent, ()))
    memo[term] = result
    return result


def propagate_scores(
    terms: Sequence[str], scores: Sequence[float], active: set[str], parents: Mapping[str, set[str]]
) -> Dict[str, float]:
    propagated: Dict[str, float] = {}
    memo: Dict[str, set[str]] = {}
    for term, raw_score in zip(terms, scores):
        if term not in active:
            continue
        score = float(raw_score)
        propagated[term] = max(score, propagated.get(term, 0.0))
        for parent in ancestors(term, parents, memo):
            if parent in active:
                propagated[parent] = max(score, propagated.get(parent, 0.0))
    return propagated


def find_structure(structure_dir: Path | None, protein_id: str) -> Path | None:
    if structure_dir is None:
        return None
    names = (
        f"AF-{protein_id}-F1-model_v6.pdb.gz", f"AF-{protein_id}-F1-model_v6.pdb",
        f"AF-{protein_id}-F1-model_v4.pdb.gz", f"AF-{protein_id}-F1-model_v4.pdb",
        f"{protein_id}.pdb.gz", f"{protein_id}.pdb",
    )
    return next((structure_dir / name for name in names if (structure_dir / name).is_file()), None)


def structure_edges(path: Path, expected_length: int) -> Tuple[np.ndarray, np.ndarray] | None:
    coords: List[Tuple[float, float, float]] = []
    first_chain = None
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[16:17] not in (" ", "A"):
                continue
            chain = line[21:22]
            if first_chain is None:
                first_chain = chain
            if chain != first_chain:
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if len(coords) != expected_length:
        return None
    xyz = np.asarray(coords, dtype=np.float32)
    delta = xyz[:, None, :] - xyz[None, :, :]
    row, col = np.where(np.einsum("ijk,ijk->ij", delta, delta) <= 100.0)
    return row.astype(np.int64), col.astype(np.int64)


def sequence_edges(contacts, length: int, threshold: float) -> Tuple[np.ndarray, np.ndarray]:
    contact_array = contacts[:length, :length].detach().float().cpu().numpy()
    row, col = np.where(contact_array >= threshold)
    if length > 1:
        backbone = np.arange(length - 1, dtype=np.int64)
        row = np.concatenate((row, backbone, backbone + 1))
        col = np.concatenate((col, backbone + 1, backbone))
    return row.astype(np.int64), col.astype(np.int64)


def build_graph(sequence: str, embedding, edges: Tuple[np.ndarray, np.ndarray]):
    import dgl
    import torch

    length = len(sequence)
    row, col = edges
    graph = dgl.graph((torch.from_numpy(row), torch.from_numpy(col)), num_nodes=length)
    graph = dgl.add_self_loop(dgl.remove_self_loop(graph))
    graph.ndata["feat"] = embedding[:length].detach().cpu()
    graph.ndata["native_x"] = torch.tensor(
        [AA_TO_INDEX.get(aa, AA_TO_INDEX["X"]) for aa in sequence], dtype=torch.int32
    )
    return graph


def make_model(aspect: str, part: str, model_dir: Path, device):
    import torch
    from megraph.layers.graph_layers.layers import GCNLayer
    from megraph.models.megraph import MeGraph

    def build_conv(**kwargs):
        return GCNLayer(**kwargs)

    part_index = {"All": 0, "Head": 1, "Tail": 2}[part]
    output_dim = OUTPUT_DIMS[aspect][part_index]
    model = MeGraph(
        input_dims=[0, 512, 0], output_dims=[output_dim, 0, 0], pe_dim=0, task="gpred",
        build_conv=build_conv, n_layers=3, g_hidden=0, n_hidden=512, e_hidden=0,
        activation="relu", dropout=0.2, norm_layer="layer", stem_beta=1, branch_beta=1,
        keep_beta=1, cross_beta=1, last_hidden_dims=[2048], max_height=1,
        lambdaw=0.2, eps=0.2, temperature=0.15, pool_aggr_edge="sum",
        cluster_size_limit=3, perturbed=False, soft_readout=True,
        global_pool_methods=["mean"],
    ).to(device)
    checkpoint_path = model_dir / f"best_{aspect}_{part}.pt"
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"Missing checkpoint: {checkpoint_path}")
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    model.set_perturbed(False)
    return model


def load_models(model_dir: Path, device):
    return {
        aspect: {part: make_model(aspect, part, model_dir, device) for part in ("All", "Head", "Tail")}
        for aspect in ASPECTS
    }


def specialized_indices(project_dir: Path, aspect: str) -> Tuple[np.ndarray, np.ndarray]:
    path = project_dir / "data" / "DMETrain" / f"distribution_{aspect}_All.txt"
    counts = np.asarray([int(value) for value in path.read_text().replace(" ", "").split(",") if value.strip()])
    expected = OUTPUT_DIMS[aspect][0]
    if len(counts) != expected:
        raise ValueError(f"Expected {expected} values in {path}, found {len(counts)}")
    head_count = OUTPUT_DIMS[aspect][1]
    ranked = np.argsort(-counts, kind="stable")
    return np.sort(ranked[:head_count]), np.sort(ranked[head_count:])


def predict_aspect(graph, models, head_indices: np.ndarray, tail_indices: np.ndarray, device) -> np.ndarray:
    import torch

    graph = graph.to(device)
    outputs = {}
    with torch.inference_mode():
        for part, model in models.items():
            # MeGraph writes its projected 512-dimensional representation back to
            # graph.ndata["feat"]. Isolate each ensemble member so the next model
            # still receives the original 1280-dimensional ESM embedding.
            with graph.local_scope():
                primary, label_graph = model(graph)
                outputs[part] = torch.sigmoid(0.5 * primary + 0.5 * label_graph)
        combined = outputs["All"].clone()
        combined[:, head_indices] = 0.5 * (combined[:, head_indices] + outputs["Head"])
        combined[:, tail_indices] = 0.5 * (combined[:, tail_indices] + outputs["Tail"])
    return combined[0].float().cpu().numpy()


def write_predictions(path: Path, rows: Iterable[Tuple[str, str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
            for protein_id, term, score in rows:
                handle.write(f"{protein_id}\t{term}\t{score:.6f}\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    query_path, graph_path, output_path = map(Path, (args.query_file, args.graph, args.output_file))
    for path in (query_path, graph_path):
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    import dgl
    import esm
    import torch

    torch.set_num_threads(args.num_threads)
    if args.device == "auto":
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available")

    project_dir = Path(__file__).resolve().parent
    vocabulary = load_term_vocabulary(project_dir / "data" / "nrPDB-GO_2019.06.18_annot.tsv")
    active_terms, parents = load_ontology(graph_path)
    indices = {aspect: specialized_indices(project_dir, aspect) for aspect in ASPECTS}
    models = load_models(Path(args.model_dir), device)
    esm_checkpoint = Path(args.esm_model)
    regression_checkpoint = esm_checkpoint.with_name(esm_checkpoint.stem + "-contact-regression.pt")
    if not esm_checkpoint.is_file() or not regression_checkpoint.is_file():
        raise FileNotFoundError(
            "ESM-1b weights are not mounted. Expected both "
            f"{esm_checkpoint} and {regression_checkpoint}"
        )
    esm_model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(esm_checkpoint))
    esm_model = esm_model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()
    structure_dir = Path(args.structure_dir) if args.structure_dir else None

    record_iterator = iter(read_fasta(query_path))
    first_record = next(record_iterator, None)
    if first_record is None:
        raise ValueError(f"No FASTA records found in {query_path}")
    records = chain((first_record,), record_iterator)
    seen: set[str] = set()

    def prediction_rows():
        for number, (protein_id, raw_sequence) in enumerate(records, 1):
            if protein_id in seen:
                raise ValueError(f"Duplicate query identifier: {protein_id}")
            seen.add(protein_id)
            sequence = "".join(aa if aa in AA_TO_INDEX else "X" for aa in raw_sequence)[: args.max_length]
            if not sequence:
                raise ValueError(f"Empty sequence for {protein_id}")
            _, _, tokens = batch_converter([(protein_id, sequence)])
            with torch.inference_mode():
                result = esm_model(tokens.to(device), repr_layers=[33], return_contacts=True)
            embedding = result["representations"][33][0, 1 : len(sequence) + 1]
            structure = find_structure(structure_dir, protein_id)
            edges = structure_edges(structure, len(sequence)) if structure else None
            edge_source = "structure" if edges is not None else "ESM contacts"
            if edges is None:
                edges = sequence_edges(result["contacts"][0], len(sequence), args.contact_threshold)
            graph = build_graph(sequence, embedding, edges)
            all_scores: Dict[str, float] = {}
            for aspect in ASPECTS:
                scores = predict_aspect(graph, models[aspect], *indices[aspect], device)
                for term, score in propagate_scores(vocabulary[aspect], scores, active_terms, parents).items():
                    if term != ROOT_TERMS[aspect] and score >= args.min_score:
                        all_scores[term] = max(score, all_scores.get(term, 0.0))
            ranked = sorted(all_scores.items(), key=lambda item: (-item[1], item[0]))
            if args.n_terms:
                ranked = ranked[: args.n_terms]
            for term, score in ranked:
                yield protein_id, term, min(1.0, max(0.0, score))
            print(f"[{number}] {protein_id}: {len(ranked)} terms ({edge_source})", file=sys.stderr)

    write_predictions(output_path, prediction_rows())
    print(f"Predictions written to {output_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"GOBoost failed: {error}", file=sys.stderr)
        raise
