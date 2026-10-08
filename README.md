# Coherence Response Index (CRI)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23229118.svg)](https://doi.org/10.5281/zenodo.23229118)

Code for the article

> Hidenori Tani. **The Coherence Response Index: predicted microRNA targets show no coherent regulation in chemically exposed mouse embryonic stem cells.** *Chemico-Biological Interactions* (2026). https://doi.org/10.1016/j.cbi.2026.112371

The CRI asks a simple question of a miRNA expression change: do the miRNA's *predicted* targets actually move in the expected direction in matched mRNA data from the same samples? For each responsive miRNA it compares the fold changes of its predicted targets with an expression-matched background of non-target genes (Mann–Whitney AUC, Benjamini–Hochberg FDR) and classifies the response as productive, anti-coherent or passenger. This repository contains the Python package that implements the test, the scripts that re-quantify the public data and run every analysis in the article, and the scripts that draw the figures and Table 1.

## Contents

| Path | What it is |
|---|---|
| `src/tmc/` | The package (`tmc`): accession map (`config`), expression input and the two fold-change estimators (`io_expr`), target tables (`targets`), expression-matched background (`nullmodel`), coherence test and CRI (`coherence`), generic vs chemical-specific structure (`specificity`), support of predicted-target pathway claims by observed repression (`audit`), end-to-end orchestration (`pipeline`) |
| `tests/` | 33 unit tests (`pytest`) |
| `scripts/` | Data retrieval, quantification and every analysis (see the order below) |
| `figures/` | Figures 1–5 (`make_figs_R1.py`) and Table 1 (`make_tables.py`); `figstyle.py` is the shared plotting style |
| `config/` | Run accessions with ENA md5 (`accessions.tsv`, `posctrl_accessions.tsv`), published anchors used in the reproducibility check (`published_top30.tsv`), RefSeq → gene symbol map from GENCODE vM36 (`refseq2symbol.tsv`) |
| `data/*/SOURCE.txt`, `data/*/SHA256SUMS.txt` | Where and when the two manually downloaded inputs were obtained, with their checksums |

## Requirements

- Python ≥ 3.11 with numpy, scipy, pandas, matplotlib, statsmodels, pyyaml (installed by `pip install -e .`); `pytest` for the tests. Released and tested with Python 3.13.5, numpy 2.4.4, scipy 1.17.1, pandas 2.3.3, matplotlib 3.10.8, statsmodels 0.14.6.
- For re-quantification only: `kallisto`, `bowtie2`, `cutadapt` and `curl` on the `PATH` (versions on the analysis machine at release: kallisto 0.52.0, bowtie2 2.5.5, cutadapt 5.2).
- Disk: the long-RNA layer is about 48 GB of FASTQ in total; the scripts stream one run at a time and stop if free space falls below 6 GB.

```bash
pip install -e .
python -m pytest            # 33 tests
```

## Data

All sequencing data are public and are downloaded by the scripts (FASTQ from ENA, checked against the md5 in `config/`):

- Paired toxicogenomic data from the same RNA preparations of mouse embryonic stem cells exposed for 24 h to nine chemicals: small-RNA-seq DDBJ/ENA DRR127916–DRR127935 and long-RNA-seq DRR082819–DRR082838.
- miRNA-mimic transfection positive control: GEO GSE269995.
- Dicer1-knockout positive control: GEO GSE139348 (gene-level counts deposited by the authors).

Target predictions are **not redistributed here**; the scripts rebuild them from the official releases:

- miRDB v6.0 (`scripts/build_targets.py` downloads it).
- TargetScan Mouse 8.0: download `Predicted_Targets_Info.default_predictions.txt.zip` and `miR_Family_Info.txt.zip` from https://www.targetscan.org/mmu_80/mmu_80_data_download/ into `data/targetscan_mmu80/`, check them against `data/targetscan_mmu80/SHA256SUMS.txt`, and unzip them there. The rebuilt `config/targetscan_mmu.tsv` should have SHA-256 `5b2c6adfb08c72e85058bd0e9ca33e67b4b9f6619f1483d87d7230a3c6d5ee29`.
- For the Dicer1 control, download `GSE139348_mESCs_dicer_RNAseq_counts_DEseq.csv.gz` from the GEO supplementary folder named in `data/physctrl_GSE139348/SOURCE.txt` into that directory and check it against `SHA256SUMS.txt`.

References: GENCODE vM36 mouse transcriptome and metadata, and miRBase mature sequences (the current `mature.fa`), fetched by `scripts/fetch_refs.py` into `ref/`.

## Reproducing the analyses

Run from the repository root, in this order. Outputs go to `results/` (and `results/revision_r1/`).

```bash
python scripts/fetch_refs.py            # GENCODE vM36 kallisto index, transcript-to-gene map, miRBase bowtie2 index
python scripts/build_targets.py         # miRDB v6.0 -> config/mirdb_mmu.tsv
python scripts/build_targetscan.py      # TargetScan Mouse 8.0 -> config/targetscan_mmu*.tsv
python scripts/quant_mirna.py           # small RNA: cutadapt + bowtie2 -> results/mirna_counts.tsv
python scripts/quant_mrna.py            # long RNA: kallisto -> results/gene_tpm.tsv
python scripts/build_log2fc.py          # per-chemical fold changes (robust pseudocount and zero-replacement estimators)
python scripts/run_all.py all           # primary coherence test (results/coherence.tsv) and specificity
python scripts/sanity_vs_published.py   # recomputed vs published biomarker fold changes
python scripts/floor_demo.py            # zero-replacement floor artifact on a single quantification
python scripts/posctrl_quant.py         # mimic-transfection control: quantification
python scripts/posctrl_cri.py           # mimic-transfection control: CRI
python scripts/physctrl.py              # Dicer1-knockout control
python scripts/sensitivity.py           # robustness variants
python scripts/power_check.py           # spike-in power analysis
python scripts/revision_r1.py all       # replicate halves, background windows, CRI weighting, TargetScan, gene-type backgrounds and other variants
python figures/make_tables.py           # Table 1
python figures/make_figs_R1.py          # Figures 1-5 (written to submission_R1/fig/)
```

The background draws use Python's `random` module with seed 0, so the results are deterministic; `revision_r1.py` first re-runs the primary analysis in memory and checks that it reproduces `results/coherence.tsv` before running any variant. The derived numbers reported in the article are also in its Supplementary Tables S1–S3.

`figures/figstyle.py` is a copy of the author's general plotting-style helper with comments and docstrings removed; with it, `make_figs_R1.py` reproduces the figure PNGs submitted with the accepted manuscript pixel for pixel (checked at release).

## Use of AI assistance

As declared in the article, the analysis code (quantification wrappers, coherence statistics and matched null model, robustness, replicate and power analyses, figure scripts and unit tests) was written with Claude Code (Anthropic) to the author's specifications. The author designed the CRI and fixed its decision rule before the real-data run, set every analytical threshold, and verified the code through the unit tests, deterministic re-runs and inspection of intermediate results.

## Citation

Please cite the article above. The code is archived on Zenodo: https://doi.org/10.5281/zenodo.23229118 (all versions; v1.0.0, the code used for the article, is https://doi.org/10.5281/zenodo.23229119). `CITATION.cff` gives the citation in machine-readable form.

## License

MIT (see `LICENSE`).
