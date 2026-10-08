# Data provenance and reuse check

The supplied analysis inputs were copied byte-for-byte from the author's
earlier local reproducibility package. Their SHA-256 digests are recorded in
`null_sensitivity.json` and checked by `code/verify_release.py`. This release
does not include the scripts or a documented chain that generated the binary
KO matrix from original genome annotations. Accordingly, it makes only the
**final permutation calibration** directly rerunnable, not the entire
upstream analysis.

- **HRGM2 catalogue:** Ma, J., Kim, N., Cha, J.H. et al. *A human gut
  metagenome-assembled genome catalogue spanning 41 countries supports
  genome-scale metabolic models*. Nature Microbiology 11, 317–334 (2026).
  DOI: https://doi.org/10.1038/s41564-025-02206-1. The HRGM2 download site
  is https://www.decodebiome.org/HRGM2/. Check the source dataset's attribution
  and applicable licence before redistributing derived records.
- **GTDB taxonomy:** Release r220; https://gtdb.ecogenomic.org/. GTDB states
  that its database is offered under CC BY-SA 4.0:
  https://gtdb.ecogenomic.org/licenses. Attribute GTDB and consider whether
  share-alike applies to any modified taxonomy table being distributed.
- **KEGG Orthology:** KO identifiers and textual annotations in the matrix
  and observed association table derive from KEGG-based annotation. KEGG's
  terms are at https://www.kegg.jp/kegg/legal.html. We have **not** established
  that those terms permit unrestricted public redistribution of the derived
  matrix or the annotation-rich table. Author/institutional review is needed
  before these files are posted publicly.

Files and roles:

| File | Role | Upstream derivation documented here? |
|---|---|---|
| `data/hyptest_X.pkl` | 693-species binary KO analysis matrix | No |
| `data/hrgm2_gtdb_r220.tsv` | HRGM2-to-GTDB r220 taxon labels | No |
| `data/amr_ko_set.csv` | AMR KO set for richness adjustment | No |
| `output/tables/relabelled_associations.csv` | Observed pairs and annotation/curation labels | No |

No licence is asserted here for original code, derived tables, or external
databases. The author should choose a code licence separately if desired and
should not infer data-reuse permission from a code licence.
