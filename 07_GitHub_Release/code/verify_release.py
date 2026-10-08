"""Fast checks of the published 100-run result and the exact input bytes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from statistics import mean, stdev

ROOT = Path(__file__).resolve().parent.parent
result = json.loads((ROOT / "null_sensitivity.json").read_text(encoding="utf-8"))

for relative_path, expected in result["source_sha256"].items():
    path = ROOT / relative_path
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    assert digest == expected, f"SHA-256 mismatch: {relative_path}"

assert result["replicates"] == len(result["runs"]) == 100
assert result["n_species"] == 693
assert result["n_genes"] == 5313
assert result["n_all_pairs"] == 14_111_328
assert sum(result["pair_counts_by_support"].values()) == result["n_all_pairs"]
assert sum(result["pair_counts_by_support"][band] for band in ("4", "5-9", "10+")) == 2_394_412
assert sum(result["observed_screen_counts"][band] for band in ("4", "5-9", "10+")) == 3_034

for label in ("fixed_richness", "recomputed_richness"):
    for band in (*result["pair_counts_by_support"], "supported_4_plus"):
        counts = [run[label][band] for run in result["runs"]]
        summary = result[label + "_summary"][band]
        assert counts == summary["permutation_counts"]
        assert abs(mean(counts) - summary["mean_false_pairs"]) < 1e-10
    for run in result["runs"]:
        assert run[label]["supported_4_plus"] == sum(run[label][b] for b in ("4", "5-9", "10+"))

exact_four = result["recomputed_richness_summary"]["4"]["mean_false_pairs"]
pooled = result["recomputed_richness_summary"]["supported_4_plus"]["mean_false_pairs"]
assert exact_four == 13.25 and pooled == 18.94
assert result["fixed_richness_summary"]["supported_4_plus"]["mean_false_pairs"] == 18.96
sd = stdev(run["recomputed_richness"]["supported_4_plus"] for run in result["runs"])
assert round(sd, 2) == 3.63

print("PASS: four exact input hashes and all 100 recorded null runs")
print(f"Exactly four genera: {exact_four}/283 = {100 * exact_four / 283:.2f}%")
print(f"Four or more genera: {pooled}/3034 = {100 * pooled / 3034:.2f}% (SD {sd:.2f})")
