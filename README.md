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

**Step 1 — Inference (no HLA used as input at run time).**
HLAGuessr and THNet were run independently on the ENA/MiXCR-processed repertoire, using only TCR sequence features. This was performed before HLA ground truth was obtained for this project.

**Important caveat on HLAGuessr**: Rosati et al.'s cohort was one of three datasets (alongside Russell et al. and Emerson et al.) used in the original training of HLAGuessr's classifier (Ruiz Ortega et al. 2025), with HLA genotypes obtained directly from the authors for that purpose. This means HLAGuessr's predictions on this cohort are **not fully blind** — the model had prior exposure to (a subset of) these patients' true HLA during training. High accuracy here was therefore *expected by design*, and this step served to confirm that expectation once ground truth became available to us. **THNet, by contrast, was trained on an entirely independent set of 4,144 donors (Pan et al. 2025) and never saw this cohort** — its validation below is genuinely blind and is the more informative result for assessing real-world inference accuracy on new repertoires.

**Step 2 — Validation against ground truth.**
Once Dr. Elisa Rosati shared the HIBAG-imputed HLA genotypes, each method's predictions were compared independently against real HLA to measure accuracy. **HLA ground truth values themselves are not published in this repository** — only aggregate validation metrics are reported, in agreement with the data sharing terms.

## Repository structure
rosati_hlaguessr.ipynb                      # MiXCR processing + HLAGuessr inference

rosati_thnet_inference.ipynb                # THNet inference + method comparison

rosati_master_table_construction.ipynb      # Unpack RDS, build ground-truth master table

rosati_qc_hla_validation.ipynb              # QC filtering + validation vs HIBAG ground truth

The two validation notebooks are shared without cell outputs, since intermediate outputs would partially expose individual patient HLA genotypes.

## Results

### Quality control

A productive-CDR3 filter (`C...F/W`) and CDR1/CDR2-availability filter were applied to the 19,820,069 raw clonotypes shared by Elisa Rosati. Only **1.02% were removed**; all 192 patients retained ≥1,300 clones per chain, well above standard depth thresholds.

### Validating inference against real HLA

Predictions from each method (independently, ≥90% confidence threshold) were compared against the HIBAG genotypes, restricted to HLAGuessr's 94 modelable alleles:

| Metric | HLAGuessr ≥90% (semi-informed*) | THNet ≥90% (genuinely blind) |
|---|---|---|
| Patients validated | 175 | 178 |
| Predictions made | 1,214 | 1,024 |
| True positives | 1,188 | 1,009 |
| False positives | 26 | 15 |
| Precision (PPV) | 97.9% | 98.5% |
| Sensitivity | 56.3% | 47.1% |

*\*Rosati was part of HLAGuessr's training data — see caveat above.*

**Key finding**: THNet's result is the more meaningful benchmark for real-world performance, since it had no prior exposure to this cohort. **47.1% sensitivity at 98.5% precision on a genuinely unseen cohort is a strong result** for blind TCR-based inference, where the model relies entirely on public-clonotype signal learned from other datasets. HLAGuessr's higher sensitivity (56.3%) is consistent with — and partly explained by — its prior training exposure to this cohort, so the comparison between the two numbers should not be read as "HLAGuessr is simply the better method."

Separately, we tested requiring agreement between both methods (consensus ≥90% on both) and found this **substantially worse** than either method alone (sensitivity 8.1%, vs. 56.3%/47.1% individually) — confirming that requiring cross-method agreement discards valid signal without meaningfully improving the already-high precision (~98%) of either method on its own.

Neither result is comparable to >90% benchmarks from direct molecular HLA typing (NGS/SSO), which sequences genomic DNA rather than inferring genotype from immune repertoire signal.

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

- Ruiz Ortega et al. 2025 — HLAGuessr: *PLOS Computational Biology*, DOI: 10.1371/journal.pcbi.1012724
- Pan et al. 2025 — THNet: github.com/Mia-yao/THNet
- Rosati et al. 2022 — dataset: *Gut*, DOI: 10.1136/gutjnl-2021-325373, ENA: PRJEB50045

## Acknowledgements

HLA ground truth data kindly shared by Dr. Elisa Rosati.
