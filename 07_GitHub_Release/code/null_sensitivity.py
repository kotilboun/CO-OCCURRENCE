"""Whole-universe permutation sensitivity for fixed versus recomputed richness.

Reads the release inputs without changing them. Every permutation
independently shuffles each retained gene within each genus, then evaluates all
14,111,328 unordered screened gene pairs (rather than sampling pairs). The two
analyses use the same permuted matrix but either the observed species richness
vector or one recomputed from that matrix. Run with --replicates N; output is JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg.blas import sgemm

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent
DATA = SOURCE / "data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--replicates", type=int, default=100)
parser.add_argument("--seed", type=int, default=20260930)
parser.add_argument("--output", type=Path, default=SOURCE / "rerun_null_sensitivity.json")
args = parser.parse_args()

X = pd.read_pickle(DATA / "hyptest_X.pkl")
X = X[~X.index.duplicated(keep="first")]
tax = pd.read_csv(DATA / "hrgm2_gtdb_r220.tsv", sep="\t").set_index("HRGMv2").reindex(X.index)
genus = tax["Genus"].str.replace("g__", "", regex=False).to_numpy()
codes = pd.factorize(genus)[0]
groups = [np.flatnonzero(codes == k) for k in range(codes.max() + 1)]
n = len(X)
prev = X.sum(axis=0)
kept = prev[(prev >= 10) & (prev <= n - 10)].index.tolist()
seen, reps = set(), []
for ko in kept:
    signature = X[ko].to_numpy().tobytes()
    if signature not in seen:
        seen.add(signature)
        reps.append(ko)
multi = pd.Series(genus).value_counts()
multi = multi[multi >= 2].index
variation = (X[reps].groupby(genus).max().loc[multi] - X[reps].groupby(genus).min().loc[multi]) > 0
genes = [ko for ko in reps if variation[ko].any()]
gene_index = np.array([kept.index(ko) for ko in genes])
amr = set(pd.read_csv(DATA / "amr_ko_set.csv").iloc[:, 0])
rich_index = np.array([j for j, ko in enumerate(kept) if ko not in amr])
M = X[kept].to_numpy(dtype=np.float32)
V = np.asfortranarray(variation[genes].to_numpy(dtype=np.float32))
G = len(genes)
upper = np.triu_indices(G, 1)
support = sgemm(1.0, V, V, trans_a=True)[upper].astype(np.int16)
bands = [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 3, "3"),
         (4, 4, "4"), (5, 9, "5-9"), (10, 1000, "10+")]
masks = {label: (support >= lo) & (support <= hi) for lo, hi, label in bands}
pair_counts = {label: int(mask.sum()) for label, mask in masks.items()}


def centered(A: np.ndarray) -> np.ndarray:
    result = A.copy()
    for ix in groups:
        result[ix] -= A[ix].mean(axis=0)
    return result


def correlations(A: np.ndarray, richness: np.ndarray) -> np.ndarray:
    C = centered(A)
    r = centered(richness[:, None])[:, 0]
    ss = r @ r
    if ss > 0:
        C -= np.outer(r, (r @ C) / ss)
    scale = C.std(axis=0)
    scale[scale == 0] = np.inf
    C /= scale * np.sqrt(n)
    C = np.asfortranarray(C)
    return sgemm(1.0, C, C, trans_a=True)[upper]


observed_richness = M[:, rich_index].sum(axis=1)
observed_corr = correlations(M[:, gene_index], observed_richness)
observed_screen_counts = {
    band: int((np.abs(observed_corr[mask]) >= 0.5).sum())
    for band, mask in masks.items()
}
association_table = pd.read_csv(SOURCE / "output" / "tables" / "relabelled_associations.csv")
reported_screen_counts = {
    label: int(((association_table.genus_support >= lo) &
                (association_table.genus_support <= hi)).sum())
    for lo, hi, label in bands
}
print("Observed full-screen counts:", observed_screen_counts, flush=True)
print("Reported association counts:", reported_screen_counts, flush=True)
assert observed_screen_counts == reported_screen_counts
rng = np.random.default_rng(args.seed)
records = []
for rep in range(args.replicates):
    P = M.copy()
    for ix in groups:
        if len(ix) > 1:
            order = np.argsort(rng.random((len(ix), len(kept))), axis=0)
            P[ix] = M[ix[order], np.arange(len(kept))]
            assert np.array_equal(P[ix].sum(axis=0), M[ix].sum(axis=0))
    permuted_richness = P[:, rich_index].sum(axis=1)
    record = {"replicate": rep + 1}
    for label, r in (("fixed_richness", observed_richness),
                     ("recomputed_richness", permuted_richness)):
        corr = correlations(P[:, gene_index], r)
        record[label] = {
            band: int((np.abs(corr[mask]) >= 0.5).sum())
            for band, mask in masks.items()
        }
        record[label + "_max_abs_r"] = {
            band: float(np.abs(corr[mask]).max())
            for band, mask in masks.items()
        }
        record[label + "_hits_0_9"] = {
            band: int((np.abs(corr[mask]) >= 0.9).sum())
            for band, mask in masks.items()
        }
        record[label]["supported_4_plus"] = sum(record[label][band] for band in ("4", "5-9", "10+"))
    records.append(record)
    print(f"Completed whole-universe permutation {rep + 1}/{args.replicates}", flush=True)

result = {
    "design": "independent within-genus shuffles of each kept gene; all unordered screened pairs evaluated",
    "source_directory": ".",
    "source_sha256": {
        path.relative_to(SOURCE).as_posix(): sha256(path)
        for path in (
            DATA / "hyptest_X.pkl",
            DATA / "hrgm2_gtdb_r220.tsv",
            DATA / "amr_ko_set.csv",
            SOURCE / "output" / "tables" / "relabelled_associations.csv",
        )
    },
    "seed": args.seed,
    "replicates": args.replicates,
    "n_species": n,
    "n_genes": G,
    "n_all_pairs": G * (G - 1) // 2,
    "pair_counts_by_support": pair_counts,
    "observed_screen_counts": observed_screen_counts,
    "runs": records,
}
for label in ("fixed_richness", "recomputed_richness"):
    result[label + "_summary"] = {
        band: {
            "mean_false_pairs": float(np.mean([row[label][band] for row in records])),
            "permutation_counts": [row[label][band] for row in records],
        }
        for band in list(pair_counts) + ["supported_4_plus"]
    }
out = args.output
out.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(f"Wrote {out}", flush=True)
