# HLA Inference from TCRαβ Repertoires — Rosati 2022 Cohort

Validation of two independent HLA inference methods (HLAGuessr, THNet) against real HLA genotypes, using the Rosati et al. 2022 Crohn's Disease cohort.

## Background

[Rosati et al. 2022, *Gut*](https://doi.org/10.1136/gutjnl-2021-325373) (ENA: PRJEB50045) published bulk TCRαβ repertoires for four clinical groups: Crohn's Disease (CD), Ulcerative Colitis (UC), Colorectal Cancer (CRC), and Healthy controls.

**Note on nomenclature**: "CD" refers to Crohn's Disease, not Celiac Disease.

The publicly available ENA data does not include HLA genotypes. We contacted the corresponding author, **Dr. Elisa Rosati**, who shared:

- TCRα and TCRβ repertoires (RDS format) for three groups: CD, Healthy, and UC (labelled "colitis" in the source files)
- HIBAG-imputed HLA genotypes (SNP-array based) — **only for CD and Healthy patients**

This made it possible to validate HLA inference against real ground truth for the first time in this project, rather than relying solely on agreement between independent prediction methods.

## What we found in the shared data

- The shared RDS files store TRA and TRB as **separate per-patient tables**, each named with the patient ID (e.g. `1.CD.3.Blood.bulk`). A naive unpacking step overwrites the alpha-chain file with the beta-chain file if both chains share the same output filename — this was identified and corrected by writing chain-specific filenames (`_TRA.tsv` / `_TRB.tsv`).
- HLA genotypes were provided only for the CD and Healthy groups (192 patients total: 98 CD + 94 Healthy). The UC ("colitis") group has TCR data but no HLA ground truth.
- We confirmed the shared TCR data corresponds to the same sequencing runs as the publicly available ENA dataset (PRJEB50045) by cross-referencing patient IDs and verifying clonotype overlap (16–20/20 top clonotypes matched per patient checked).

## Approach

**Step 1 — Blind inference (no HLA used as input).**
HLAGuessr and THNet were run independently on the ENA/MiXCR-processed repertoire exactly as for any cohort without known HLA — using only TCR sequence features. This blind inference was performed before HLA ground truth was available.

**Step 2 — Validation against ground truth.**
Once Dr. Elisa Rosati shared the HIBAG-imputed HLA genotypes, each method's blind predictions were compared independently against real HLA to measure accuracy. **HLA ground truth values themselves are not published in this repository** — only aggregate validation metrics are reported, in agreement with the data sharing terms.

## Repository structure
01_rosati_master_table_construction.ipynb   # Unpack RDS, build ground-truth master table

02_rosati_qc_hla_validation.ipynb           # QC filtering + validation vs HIBAG ground truth

Both notebooks are shared without cell outputs, since intermediate outputs would partially expose individual patient HLA genotypes.

## Results

### Quality control

A productive-CDR3 filter (`C...F/W`) and CDR1/CDR2-availability filter were applied to the 19,820,069 raw clonotypes shared by Elisa Rosati. Only **1.02% were removed**; all 192 patients retained ≥1,300 clones per chain, well above standard depth thresholds.

### Validating blind inference against real HLA

Predictions from each method (independently, ≥90% confidence threshold) were compared against the HIBAG genotypes, restricted to HLAGuessr's 94 modelable alleles:

| Metric | HLAGuessr ≥90% | THNet ≥90% |
|---|---|---|
| Patients validated | 175 | 178 |
| Predictions made | 1,214 | 1,024 |
| True positives | 1,188 | 1,009 |
| False positives | 26 | 15 |
| Precision (PPV) | 97.9% | 98.5% |
| Sensitivity | **56.3%** | 47.1% |

**Key finding**: both methods independently achieve strong, comparable precision (~98%) when confident. HLAGuessr detects a somewhat larger fraction of true HLA alleles (56.3% vs 47.1% sensitivity), while THNet is marginally more precise. Both substantially outperform a strict dual-method consensus requirement (which we tested separately and found overly conservative — see `02_rosati_qc_hla_validation.ipynb`), confirming that requiring agreement between independent methods discards valid signal without meaningfully improving precision.

56% sensitivity at 98% precision is a strong result for *blind* TCR-based inference, where the model never observes HLA labels and must rely entirely on public-clonotype signal — this is not directly comparable to >90% benchmarks from direct molecular HLA typing (NGS/SSO), which sequences genomic DNA rather than inferring genotype from immune repertoire signal.

### Known artefacts

- DPA1*01:03 (THNet): predicted positive in nearly 100% of patients regardless of repertoire — excluded from all validation, consistent with class imbalance in THNet's training data.

## Requirements

```bash
conda create -n bioinf python=3.12
conda activate bioinf
pip install hlaguessr thnet pandas numpy scikit-learn openpyxl
```

## Citation

If you use this pipeline, please cite:

- Ortega et al. 2025 — HLAGuessr: *PLOS Computational Biology*, DOI: 10.1371/journal.pcbi.1012724
- Pan et al. 2025 — THNet: github.com/Mia-yao/THNet
- Rosati et al. 2022 — dataset: *Gut*, DOI: 10.1136/gutjnl-2021-325373, ENA: PRJEB50045

## Acknowledgements

HLA ground truth data kindly shared by Dr. Elisa Rosati.
