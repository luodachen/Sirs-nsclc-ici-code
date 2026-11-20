# Analysis code -- smoking-related immune remodeling index (SIRS) in NSCLC

Supplementary analysis code for the manuscript:

> **Smoking-related immune remodeling and response to anti-PD-L1 therapy: a signature built in
> NSCLC and tested in urothelial cancer**

This package contains the complete R pipeline that produced every result, figure and supplementary
item in the manuscript, together with the dependency-installation scripts and the patched copy of
one third-party package that is required to install it on hg38.

---

## 1. Contents

```
.
|- README.md                       this file
|- R/                              the analysis pipeline (see section 4 for the two run orders)
|   01_Step1_exposure_quantification.R        smoking exposure quantification (Step 1) - minimal set
|   02_Step1b_strategyA_CSI_robustness.R      exposure QC and parameter robustness - supplies Table 1 baseline values (supplementary tables)
|   03_Step1c_missingness_diagnosis_and_tsc.R missingness diagnosis and tsc conventions - optional
|   04_Step2_immune_microenvironment_TMB_SBS4.R  deconvolution, TMB and SBS4 (Step 2) - minimal set
|   05_Step2b_diagnostics_and_robustness.R    Step 2 diagnostics (Q1-Q4) - supplementary only
|   06_Step3_SIRS_construction.R              SIRS construction (Step 3) - minimal set
|   07_Step4_IMvigor210_external_validation.R external validation in IMvigor210 (Step 4) - minimal set
|   08_Step4b_sensitivity_analyses.R          IC0-stratum robustness analyses - supplementary only
|   09_Step4c_single_algorithm_sensitivity.R  single-deconvolution-algorithm sensitivity - supplementary only
|   10_Step4d_TCGA_cross_cancer_replication.R internal cross-cancer replication (BLCA/HNSC/ESCA) - supplementary only
|- install/
|   install_deps.R                  installs all Step 2 dependencies (idempotent)
|   install_deconstructsigs.R       installs an hg38-capable deconstructSigs
|   smoke_test_step2.R              minimal verification of the two heavy dependency chains
|- vendor/
|   deconstructSigs-1.9.0-hg38patched/  patched package source (see its PATCH_NOTE.md)
|- docs/
|   06_data_source_audit_smoking_fields.md      why the precomputed cigarettes-per-day field was not used
```

Note: `R/` is a flat directory because the scripts source nothing from one another; each script is
self-contained and communicates through files on disk. The document under `docs/` is cited by name
in the script headers.

## 2. Requirements

- **R >= 4.3** (developed and run on R 4.5.1 with Bioconductor 3.21).
- No Rtools is required: all Bioconductor dependencies are installed as binaries, and the only
  source-built packages are pure R.
- Disk: about 1 GB for the clinical, expression and mutation data (about 3 GB more if methylation
  data are added, which this pipeline does not use).
- CRAN packages: `tidyverse`, `data.table`, `mice`, `patchwork`.

## 3. Data sources

No raw data are redistributed here. All inputs are public:

| Data | Source | How it is obtained |
|---|---|---|
| TCGA-LUAD / TCGA-LUSC clinical, expression (STAR FPKM-UQ), somatic mutation | UCSC Xena GDC hub | downloaded automatically by `R/01_...R` on first run (about 0.5-1 GB; resumable) |
| TCGA-LUAD/LUSC legacy clinical matrix | UCSC Xena | same as above |
| GDC clinical API (used only for the data-source audit) | https://api.gdc.cancer.gov | queried by section 6.0 of `R/01_...R` |
| IMvigor210 (urothelial carcinoma, anti-PD-L1, n = 348) | `IMvigor210CoreBiologies` R package, http://research-pub.gene.com/IMvigor210CoreBiologies/ | **not a download in this pipeline** - obtained manually, see below |
| BLCA / HNSC / ESCA (cross-cancer replication) | UCSC Xena GDC hub | downloaded by `R/10_...R` |

**Obtaining IMvigor210 before running script `07`.** This is the one input that no script in this
package fetches, and no installation script installs it. Script `07` requires exactly one file,
`data/IMvigor210/cds.RData` (51,233,046 bytes), and stops with an error if it is missing:

1. Download the package archive from the GENENTECH site above (package version 1.0.0). If that host
   is unreachable, any HTTP mirror of the same file is acceptable, provided the downloaded
   `cds.RData` is 51,233,046 bytes; the copy analysed in the paper was retrieved through a public
   proxy mirror and verified to be byte-for-byte identical to the record inside the official archive.
2. The archive contains `IMvigor210CoreBiologies/data/cds.RData` together with `fmone.RData`,
   `human_gene_signatures.RData`, `data.R` (the roxygen data documentation) and
   `customFunctions.R` (the authors' helper functions).
3. Copy at least `cds.RData` to `data/IMvigor210/` under your `proj_dir`. The other files are not
   read by this pipeline and may be kept for provenance.

Do **not** attempt to `library(IMvigor210CoreBiologies)`: the package was built against the retired
`DESeq` (classic) and `Biobase`, so it will not install on a current Bioconductor. Script `07`
therefore reads `cds.RData` directly with `load()` and reaches into the `CountDataSet` through
`attr()`; loading `cds` in a fresh session prints a harmless warning that the `DESeq` namespace is
unavailable and has been replaced by `.GlobalEnv`. This is expected and does not affect the analysis.

## 4. How to run

1. Create a working directory and copy `R/`, `install/`, `vendor/` and `docs/` into it.
2. **Edit `proj_dir` in section 0 of each script** (it currently reads
   `D:/Projects/smoking-index-ici`). Every script uses the same string, so a single find-and-replace
   across `R/` is enough. Scripts also assume the working directory layout created by step 1.
3. Install the heavy dependencies **before** anything else:

   ```r
   setwd("<your proj_dir>"); source("install/install_deps.R")            # idempotent
   setwd("<your proj_dir>"); source("install/install_deconstructsigs.R") # hg38 build
   Rscript install/smoke_test_step2.R                                    # ~2-5 min verification
   ```

   `install_deps.R` installs `BSgenome`, `BSgenome.Hsapiens.UCSC.hg38` (about 731 MB - the slowest
   step), `quantiseqr`, `EPIC` and `MCPcounter`. `install_deconstructsigs.R` installs the patched
   `deconstructSigs` from `vendor/`. If the SBS4 or deconvolution dependency is missing, the
   corresponding columns are silently skipped by script 04, so **run the smoke test first**.
4. Run the pipeline. The scripts share no code: each is self-contained and they communicate only
   through files on disk, so the chain can be entered at any point whose upstream outputs already
   exist. Two levels of reproduction are possible.

   **(a) Minimal set - reproduces the main-text figures and Table 3** (Figures 2, 3, 4 and 5):

   ```
   01  ->  04  ->  06  ->  07
   ```

   **(b) Full set - additionally reproduces every Online Resource (1-8):**

   ```
   01  ->  02  ->  03  ->  04  ->  05  ->  06  ->  07  ->  08  ->  09  ->  10
   ```

   Script `01` downloads the raw data on first run. The scripts that appear in (b) but not in (a) are
   side-branches: `02` (Online Resource 1), `05` (Online Resource 2), `08` (Online Resource 6), `09`
   (Online Resource 7) and `10` (Online Resource 8). None of them feeds a later analysis (see the
   "Required by" column in section 7), so omitting any of them changes no number in the main text.
   One caveat about set (a): it reproduces Figures 2-5 and Table 3, but **not** the numeric columns
   of Table 1, which come from `step1b_CSI_tertile_baseline.csv` written by `02`. Add `02` if Table 1
   must be regenerated from source.

   Ordering constraints that must be respected:

   - **`03` before `04`.** This one is a *soft* dependency: `04` reads
     `step1_output/step1c_csi_by_sample.csv` only if the file exists (`if (file.exists(...))`).
     Without it `04` still writes `step2_immune_tme.csv`, but the time-since-cessation convention
     sensitivity columns and `step2_tsc_convention_sensitivity.csv` are skipped.
   - **`07` before `08` and `07` before `09`.** Both read `step4_output/step4_analysis_table.csv`,
     which `07` writes.
   - **`07` before `09`.** `09` also reads `data/IMvigor210/deconv_imvigor.rds`, a cache that `07`
     writes on its first run.
   - **`06` before `10`.** `10` applies the frozen `step3_output/step3_SIRS_weights.csv`.

   `10` is independent of the whole Step 4 chain: it writes into `step4_output/` but reads nothing
   from `07`, `08` or `09`, so once `06` has finished it can be run in parallel with `07`. `09`,
   by contrast, cannot be parallelised with `07` (it consumes the deconvolution cache that `07`
   creates), and `08` likewise needs `07`'s analysis table first.

   No intermediate outputs are redistributed with this package (only code), so a fully fresh
   reproduction starts at `01`; the ability to enter the chain mid-way applies when the upstream
   outputs are already present on disk.

   **Inputs consumed by each script.** The order above is a reading order, not a dependency chain:
   the scripts are not sequential steps of one another. This table lists what each script actually
   reads, so a reader can tell which earlier scripts are genuinely required.

   | Script | Reads (inputs) | Hard requirement? |
   |---|---|---|
   | `01` | raw data under `data_raw/xena`, downloaded on first run | none; entry point |
      | `02` | `step1_output/smoking_index_clinical.csv` (from `01`) | none; writes `step1b_*` tables only |
   | `03` | `step1_output/smoking_index_clinical.csv` (from `01`) | none; writes `step1c_*` tables only |
   | `04` | `step1_output/smoking_index_clinical.csv` (from `01`); `step1c_csi_by_sample.csv` (from `03`) **if present** | requires `01`; `03` is optional |
   | `05` | `step2_output/step2_immune_tme.csv` (from `04`) | requires `04` |
   | `06` | `step2_output/step2_immune_tme.csv` (from `04`) | requires `04` |
   | `07` | `data/IMvigor210/cds.RData` (**manual, see section 3**); `step3_output/step3_SIRS_weights.csv` (from `06`) **if present** | requires `06` and the IMvigor210 data |
   | `08` | `step4_output/step4_analysis_table.csv` (from `07`) | requires `07` |
   | `09` | `step2_output/step2_immune_tme.csv` (from `04`), `step3_output/step3_SIRS_all_samples.csv` (from `06`), `data/IMvigor210/deconv_imvigor.rds` and `step4_output/step4_analysis_table.csv` (both from `07`) | requires `04`, `06`, `07` |
   | `10` | `step3_output/step3_SIRS_weights.csv`, `step3_SIRS_all_samples.csv` and `step3_SIRS_exposure_validity.csv` (all from `06`); raw BLCA/HNSC/ESCA data | requires `06` only |

   Three consequences that are easy to get wrong:

   - **`02` is not a prerequisite for anything.** It reads only `step1_output/smoking_index_clinical.csv`
     and writes `step1b_*` tables; no later script reads them. `07` in particular does not depend on
     `02` - it needs `06` (the frozen weights) and the IMvigor210 data. Running `02` (or `03`) before
     `07` is harmless but changes nothing in the Step 4 results.
   - **`07` carries its weights in the script itself.** The eight weights and their training means and
     standard deviations are hard-coded, following the values reported in the manuscript.
     `step3_output/step3_SIRS_weights.csv` is read only as a cross-check (Guard 4, tolerance
     `1e-6`); if the file is absent the check is skipped and `07` still runs.
   - **`09` is the most input-hungry script** and cannot be run in parallel with `07`: besides
     `07`'s analysis table it needs the deconvolution cache `data/IMvigor210/deconv_imvigor.rds`,
     which `07` creates on its first run.

## 5. Key analytical conventions

These are the decisions that a reader needs in order to reproduce the results.

- **Cigarettes per day (CPD)** is *derived*, not measured:
  `CPD = pack-years x 20 / smoking duration`, computed in two passes (duration first). The
  precomputed `cigarettes_per_day` field in the UCSC Xena GDC hub was **not** used: it is a
  deterministic rescaling of pack-years (`== pack-years x 20 / 365`; exact in 782/782 cases) and the
  corresponding GDC API field has 0% coverage in TCGA-LUAD/LUSC. The full evidence chain is in
  `docs/06_data_source_audit_smoking_fields.md`; the check is re-run at every execution and written
  to `QC_xena_cpd_audit.csv`.
- **Consequence**: `cigarette-years == pack-years x 20` is a constructive identity, so the two are
  not independent exposure axes. Downstream analyses therefore use pack-years, log(cigarette-years)
  and the comprehensive smoking index instead.
- **Comprehensive smoking index (CSI)**: `(1 - 0.5^(dur*/tau)) * 0.5^(tsc*/tau) * ln(int+1)` with
  `tsc* = max(tsc - delta, 0)` and `dur* = max(dur + tsc - delta, 0) - tsc*`. The primary column
  `csi_lc` uses the published lung-cancer parameters tau = 26 years, delta = 0.7 years; three
  published-parameter sensitivity variants are also exported.
- **TMB** counts non-synonymous SNVs plus indels divided by 38 Mb. The Xena `effect` column uses VEP
  terms and can carry several effects joined together, so it must be matched with a regular
  expression rather than `%in%` against MAF vocabulary.
- **SBS4** is fitted per sample with `deconstructSigs` against the **COSMIC v2 30-signature set**
  (`signatures.cosmic`) on the hg38 reference, with exome-to-genome trinucleotide-bias correction,
  and is computed only for samples with at least 50 SNVs. Both the relative weight (`sbs4_weight`)
  and the absolute burden (`sbs4_abs = sbs4_weight x tmb_per_mb`) are exported because the two
  behave differently in downstream analyses.
- **CYT** follows Rooney 2015 (geometric mean of GZMA and PRF1). **TIS** follows the 18-gene Ayers
  2017 signature.
- **Deconvolution**: quanTIseq (`quantiseqr`) is the primary pathway; EPIC and MCP-counter are
  comparators. All three tables are normalised through the same `normalize_deconv()` helper.
- **SIRS** is an effect-size-weighted combination of eight immune features derived from
  confounder-adjusted partial correlations; the weights are frozen in `step3_output/step3_SIRS_weights.csv`
  and are never re-estimated in the validation steps.

## 6. Known pitfalls (these cost real debugging time)

- The Xena expression matrix uses 16-character sample barcodes including the `-01A`/`-01B` suffix;
  truncating to 15 characters creates duplicate column names and must be de-duplicated.
- The gene probe map file is named `gencode.v36.annotation.gtf.gene.probemap` (lower-case `m`); its
  `id` column carries version suffixes and must be stripped on both sides before joining.
- `deconstructSigs::whichSignatures()` is a **single-sample** interface: it must be called one
  sample at a time, with `tumor.ref` as a `data.frame`.
- Loading a Bioconductor genome package (`BSgenome.Hsapiens.UCSC.hg38`) brings `GenomicRanges` into
  the search path, whose S4 generics mask `purrr::reduce`, `union`, `intersect`, `setdiff`, `slice`
  and `rename`. Namespace-qualify these calls explicitly.
- After loading such packages, build model formulas with `reformulate()` rather than `paste()`;
  `paste()` inserts spaces and produces back-quoted names containing spaces.
- Figures and console text in this pipeline are ASCII-only on purpose, so that they render correctly
  under the simple encodings used by `pdf()` on Windows.

## 7. Which script produces which manuscript item

| Script | Manuscript item | In set (a) | Required by |
|---|---|---|---|
| `01` | Online Resource 1 (Figure S1, exposure QC) | yes | `04` |
| `02` | `step1b_CSI_tertile_baseline.csv` supplies the baseline characteristics of **Table 1**; also Online Resource 1 (exposure correlation and parameter-robustness tables) | - | nothing |
| `03` | Online Resource 1 (missingness funnel and tsc conventions) | optional | nothing (`04` consumes its output only if present) |
| `04` | Figure 2; Step 2.5 diagnostic tables (Online Resource 2) and dose-response table (Online Resource 3) | yes | `05`, `06`, `09` |
| `05` | Online Resource 2 (Step 2 diagnostic tables) | - | nothing |
| `06` | Figure 3 and Figure 4; Online Resource 4 (SIRS formula and training parameters) | yes | `07`, `09`, `10` |
| `07` | Figure 5; Table 3 (IC-stratified validation results) | yes | `08`, `09` |
| `08` | Online Resource 6 (IC0 head-to-head AUC with DeLong, Firth penalised logistic, bootstrap and permutation) | - | nothing |
| `09` | Online Resource 7 (single-algorithm sensitivity: EPIC / MCP-counter / quanTIseq) | - | nothing |
| `10` | Online Resource 8 (internal cross-cancer replication, negative) | - | nothing |

Figure 1, Table 1 and Table 2 are assembled from these outputs and are not generated by a single
script. Figure 1 and Table 2 draw on the exposure and cohort outputs of `01` and `04`; the numeric
columns of Table 1 come from `step1b_CSI_tertile_baseline.csv`, so **`02` is needed to regenerate
Table 1 even though it is not part of set (a)**.

Set (a) is the minimal reproduction set defined in section 4 (the figures and Table 3); the scripts
marked `-` are not required to reproduce those items. The last column shows the reverse relation -
which script consumes a given script's output - and is the definitive answer to "must I run X before
Y?". `02` appears in no such chain: it is a leaf.

### 7.1 Output files: manuscript-facing versus internal QC

A full run (set b) writes many more files than the manuscript displays. The table above names only
the **manuscript-facing** outputs. The remaining CSVs are internal quality-control / audit records:
they are deliberately retained so that a reviewer can trace every filtering decision, but none of
them is cited in the text or reproduced as a figure or table. They are listed here so that the
package contents can be reconciled one-to-one with the manuscript.

| Script | Manuscript-facing outputs (figure / table / OR) | Internal QC outputs (not in the manuscript) |
|---|---|---|
| `01` | `smoking_index_clinical.csv` -> feeds all later steps; exposure-QC plots -> Online Resource 1 | `download_manifest.csv`, `QC_xena_cpd_audit.csv`, `missingness_table.csv`, `exposure_nonfinite_diag.csv` |
| `02` | `step1b_CSI_tertile_baseline.csv` -> **Table 1** numeric columns; exposure-correlation and parameter-robustness tables -> Online Resource 1 | `step1b_tertile_counts.csv` (cross-checked in the text but the counts appear in Table 1) |
| `03` | missingness funnel and tsc-convention plots -> Online Resource 1; `step1c_csi_by_sample.csv` -> consumed by `04` | `step1c_missing_reason.csv`, `step1c_tsc_decision.csv`, `step1c_missing_funnel.csv` (data behind OR1) |
| `04` | `step2_immune_tme.csv` -> **Figure 2** and all of Step 3/4; dose-response table -> Online Resource 3 | `step2_qc_epic_nonconverged.csv`, `step2_tsc_convention_sensitivity.csv` |
| `05` | Q1-Q4 diagnostic tables -> Online Resource 2 (their conclusions are also stated in Results 3.3) | none beyond the four OR2 tables |
| `06` | `step3_SIRS_weights.csv` and `step3_SIRS_formula.md` -> Online Resource 4 and **Table 2**; `step3_Cindex_comparison.csv` -> **Figure 4b**; `step3_SIRS_exposure_validity.csv` -> **Figure 3b** | `step3_screen_smoking_immune.csv`, `step3_lasso_diagnostics.csv` |
| `07` | `step4_analysis_table.csv` -> **Figure 5**, **Table 3** and the base of OR6/OR7; cross-platform table -> Figure 5b | `step4_q1_crossplatform_scale.csv` (plotted in Fig 5b), `07_Step4_summary.md` |
| `08` | Firth / bootstrap / DeLong / FDR-EPV tables -> Online Resource 6 | `08_Step4b_summary.md` |
| `09` | single-algorithm AUC / complementarity tables -> Online Resource 7 | `09_Step4c_summary.md`, `step4c_index_scores_imvig.csv` |
| `10` | replication / meta-analysis / dose-response / negative-control tables -> Online Resource 8 | `10_Step4d_summary.md`, per-cohort deconv caches (`step4d_deconv_*.csv`) |

In short: every script writes at least one manuscript-facing file; the extra CSVs are audit trails.
No result cited in the manuscript is generated by a script that is absent from this package. The only
external input that no script produces is `data/IMvigor210/cds.RData` (section 3).

## 8. Third-party code

`vendor/deconstructSigs-1.9.0-hg38patched/` is a patched copy of **deconstructSigs 1.9.0**
(Rosenthal et al., 2016; https://github.com/raerose01/deconstructSigs), redistributed here under
its original licence because CRAN has archived the package and the released versions cannot be used
with hg38 coordinates. The two changes applied are documented in
`vendor/deconstructSigs-1.9.0-hg38patched/PATCH_NOTE.md`; no algorithmic behaviour was modified.

## 9. Availability

No patient-identifiable data are included. All datasets are public (see section 3). The analysis
code in this package is available at the URL given in the Data availability statement of the
manuscript.
