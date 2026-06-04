# HLA Inference from TCRαβ Repertoires

Computational pipeline for HLA genotype inference from bulk T-cell receptor sequencing data, applied to the Rosati 2022 cohort (PRJEB50045). Two independent methods are benchmarked and cross-validated.

---

## Repository structure

```
├── rosati_hla_inference.ipynb        # Full pipeline: MiXCR → HLAGuessr inference
├── thnet_rosati_hla_inference.ipynb  # THNet inference + HLAGuessr vs THNet comparison
└── README.md
```

---

## Dataset

**Rosati et al. 2022** — PRJEB50045 (ENA)

| Group | N patients | Description |
|---|---|---|
| Healthy | 99 | Healthy controls |
| CD | 114 | Celiac disease (confirmed diagnosis) |
| UC | 41 | Ulcerative colitis |
| CRC | 15 | Colorectal cancer |
| **Total** | **269** | |

Raw data: bulk paired-end TCRαβ amplicon sequencing (Illumina HiSeq 2500), MiXCR `milab-human-rna-tcr-umi-multiplex` preset, UMI-corrected.

> **Note**: No ground-truth HLA typing is publicly available for this cohort. All predictions are blind inference. Validation against serological or NGS HLA data is pending.

---

## Methods

### HLAGuessr (Ortega et al. 2025)
- **Input**: TCRα + TCRβ CDR3 sequences + V-gene family
- **Approach**: Exact matching of HLA-associated public clonotypes → L1-regularised logistic classifier per allele
- **Alleles**: 94 modelable alleles (subset with sufficient training carriers)
- **Training**: 1,039 HLA-typed donors (includes Rosati Healthy + CD subsets)
- **Mode**: `untyped=True` — blind inference, no ground truth required
- **Reference**: [doi:10.1371/journal.pcbi.1012724](https://doi.org/10.1371/journal.pcbi.1012724)

### THNet (Pan et al. 2025)
- **Input**: TCRβ CDR3 + V-gene (full IMGT notation)
- **Approach**: CDR3β + V-gene concatenation as sequence features → per-allele classifiers
- **Alleles**: 208 alleles — strict superset of HLAGuessr's 94
- **Training**: 4,144 HLA-typed donors (external to HLAGuessr training set)
- **Reference**: [github.com/Mia-yao/THNet](https://github.com/Mia-yao/THNet)

---

## Key results

### Quality control
Minimum threshold: ≥1,000 clones per chain (TRA and TRB). CRC samples fail QC (median 8 clones/patient — insufficient T-cell data). Final QC-passing cohort: **233 patients** (Healthy 96, CD 104, UC 33, CRC 3).

### HLA inference coverage

| Group | HLAGuessr (≥1 allele) | THNet (≥1 allele) | Both ≥90% |
|---|---|---|---|
| Healthy | 77/96 | 99/99 | 80/99 (80.8%) |
| CD | 83/104 | 114/114 | 82/114 (71.9%) |
| UC | 24/33 | 40/40 | 28/40 (70.0%) |
| CRC | 3/3 | 14/14 | 2/15 (13.3%) |

### Method agreement (HLAGuessr vs THNet)
Comparison restricted to 94 shared alleles (DPA1*01:03 excluded — model artefact).

| Group | Mean Jaccard similarity |
|---|---|
| Healthy | 0.526 |
| CD | 0.473 |
| UC | 0.357 |
| CRC | 0.107 |

High-confidence genotypes (both methods ≥90%): **192/267 patients**, max 14 alleles per patient (consistent with diploid HLA).

### Celiac disease — causal HLA haplotype recovery
All 114 CD patients have confirmed celiac disease. Expected causal haplotypes (DQ2.5, DQ2.2, DQ8) in >90% of cases (literature).

| | HLAGuessr | THNet | Both agree | At least one |
|---|---|---|---|---|
| DQ2.5 | 11/114 | 12/114 | 9/114 | — |
| DQ8 | 6/114 | 7/114 | 2/114 | — |
| Any causal HLA | 24/114 (21%) | 30/114 (26%) | 17/114 (15%) | 35/114 (31%) |

> **Interpretation**: The 69% gap relative to the >90% literature expectation reflects a **known methodological limitation** — bulk TCRβ repertoires have insufficient depth in class II HLA-associated public clonotypes for reliable DQA1/DQB1 inference. This is not a patient misclassification. THNet recovers 6 additional cases over HLAGuessr, suggesting higher class II sensitivity through sequence embedding vs exact clonotype matching.

### Known artefacts
- **DPA1*01:03** (THNet): predicted True in 100% of patients regardless of repertoire content. Excluded from all analyses. Likely caused by class imbalance in THNet training (~95% DPA1*01:03 frequency in European training donors).

---

## Requirements

```bash
conda create -n bioinf python=3.12
conda activate bioinf
pip install hlaguessr thnet pandas numpy scikit-learn matplotlib openpyxl
```

MiXCR v4.7+ required for raw FASTQ processing (Section 2 of `rosati_hla_inference.ipynb`).

---

## Output files

| File | Description |
|---|---|
| `rosati_hla_predictions.xlsx` | Predicted HLA genotype per patient (HLAGuessr, all loci) |
| `thnet_rosati_predictions.csv` | Raw THNet scores + calls, 267 patients × 208 alleles |
| `hlaguessr_vs_thnet_agreement.xlsx` | Allele-level agreements by group and HLA class |
| `celiac_causal_hla_hlaguessr_vs_thnet.xlsx` | Causal HLA haplotype calls for 114 CD patients |
| `rosati_hla_highconfidence_bothMethods.xlsx` | High-confidence genotype (both methods ≥90%) |

---

## Citation

If you use this pipeline, please cite:

- Ortega et al. 2025 — HLAGuessr: *PLOS Computational Biology*, DOI: 10.1371/journal.pcbi.1012724
- Pan et al. 2025 — THNet: github.com/Mia-yao/THNet
- Rosati et al. 2022 — dataset: PRJEB50045

---

## Status

- [x] MiXCR processing — 1,448 runs, 269 patients
- [x] HLAGuessr inference — 269 patients × 94 alleles
- [x] THNet inference — 267 patients × 208 alleles
- [x] Method comparison — Jaccard, agreement tables, celiac haplotype analysis
- [ ] Ground-truth HLA validation (pending response from data owners)
- [ ] Russell 2022 dataset (PRJEB52108) — next cohort
