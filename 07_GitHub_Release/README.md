# Gene co-occurrence: permutation-calibration release

This folder is the **code-and-data release for the manuscript's central
permutation-null calculation**, prepared for the author to upload to
`https://github.com/kotilboun/coocurance`. It is not the journal-submission
folder, and it does not contain the blinded manuscript, cover letter, figures,
or supplementary workbook. This release does **not** reconstruct the upstream
HRGM2/GTDB/KO annotation matrix or the literature-curation decisions from raw
genomes; those are separate stages of the study.

## What the result means

The script evaluates every one of 14,111,328 screened unordered gene pairs in
each of 100 independent within-genus permutations. There are 2,394,412
possible pairs with support from at least four genera. At `|r| >= 0.50`, the
mean null count in that pooled group is 18.94, versus 3,034 observed
associations: `18.94 / 3034 = 0.62%`. This is a **null-based point estimate**,
not a formal false-discovery-rate bound. The exactly-four-genera stratum is
different: `13.25 / 283 = 4.7%`. Pooling in higher-support strata lowers the
overall ratio; the two percentages are not contradictory.

The primary result uses richness recomputed on every permuted matrix. The
fixed-richness sensitivity has a pooled mean of 18.96. The published result
and all 100 per-permutation counts are in `null_sensitivity.json`.

## Contents

| Path | Purpose |
|---|---|
| `code/null_sensitivity.py` | Complete-pair, within-genus null calculation |
| `code/verify_release.py` | Fast integrity and numerical-consistency checks |
| `data/hyptest_X.pkl` | Analysis-ready species-by-KO binary matrix (Python pickle) |
| `data/hrgm2_gtdb_r220.tsv` | HRGM2 IDs and GTDB r220 taxonomy used in the analysis |
| `data/amr_ko_set.csv` | KO identifiers excluded from the richness covariate |
| `output/tables/relabelled_associations.csv` | Observed association table used to cross-check screen counts |
| `null_sensitivity.json` | Recorded 100-run null output |
| `DATA_PROVENANCE.md` | Source attribution and redistribution cautions |

## Check and rerun

Use Python 3.13 with the versions in `requirements.txt`. From this folder:

```bash
python -m pip install -r requirements.txt
python code/verify_release.py
python code/null_sensitivity.py --replicates 100
```

The last command is computationally substantial: it performs 100 full screens
of 14.1 million pairs. It writes `rerun_null_sensitivity.json` so that the
recorded result is not overwritten. For a quick execution smoke test, use
`--replicates 1`; that tests the pipeline but cannot estimate the 100-run
mean. The seed is 20260930. Output should match the recorded result with the
pinned environment; floating-point details may vary on other BLAS/platform
combinations. The script asserts that it recovers the observed association
counts before running permutations.

`hyptest_X.pkl` is a Python pickle and should only be opened in a trusted
environment. The release contains no raw genome sequences. This folder is a
minimal release of the final null analysis, not a claim that every manuscript
figure, curation verdict, or upstream preprocessing step is independently
rebuildable from it.

## Before a public upload

Read `DATA_PROVENANCE.md` and confirm that the author's rights and any
upstream database terms permit redistribution of the KO-derived matrix and
annotation table. **Do not publish these two files until that check is
complete.** If redistribution is not permitted, replace the files with a
rights-compliant access route and revise the README and availability statement
accordingly. Do not upload the journal manuscript to GitHub during review.

After uploading, record an immutable commit or release URL and use that URL
in the manuscript's Data and Code Availability statement. The local folder
alone does not establish that the repository is public or complete.
