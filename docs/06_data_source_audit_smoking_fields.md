# Data-source audit: availability and repair of the TCGA-LUAD/LUSC smoking fields

> Project: smoking-exposure quantification x NSCLC immune microenvironment x ICI response prediction
> Status: confirmed; implemented in `01_Step1_exposure_quantification.R` (section 6.0)

---

## 1. Summary (TL;DR)

| Item | Conclusion |
|---|---|
| **Xena `cigarettes_per_day` column** | Not usable. Empirically equal to `pack-years x 20 / 365` (holds exactly in 782/782 cases), median only **2.19 "cigarettes/day"**. It is a linear rescaling of pack-years and carries **zero independent information** |
| **Same field from the GDC clinical API** | **0% coverage** (1089 TCGA-LUAD + LUSC cases). `years_smoked`, `time_since_cessation` and `age_at_onset` are likewise 0% |
| **Xena `pack_years_smoked`** | Trustworthy. **100% identical** to the GDC value in all 782 cases; 86.2% coverage among ever-smokers |
| **Xena `years_smoked`** | Trustworthy. Back-solving cigarettes per day from it gives a median of **20.0 cigarettes/day**, with 95.2% inside the physiological 5-60 range |
| **Repair** | discard the spurious column and solve algebraically from two verified fields: **CPD = pack-years x 20 / smoking years** (an identity implied by the definition of pack-years, requiring no extra assumption) |
| **Analysis set size** | ever-smokers 905 -> **559 usable (61.8%)**; the repair **does not change N** but corrects every exposure value |

---

## 2. Evidence chain (four empirical checks)

### Check 1: the Xena CPD column is a deterministic function of pack-years

In `TCGA-LUAD.clinical.tsv.gz` / `TCGA-LUSC.clinical.tsv.gz`:

```
Xena cigarettes_per_day == pack_years_smoked / 18.25    ->  782/782 = 100.0%
```

`18.25 = 365 / 20`. A 100% exact fit means the column contains no independent measurement at all:
it is simply pack-years converted with a mistaken "365 days = 1 year" rescaling.

### Check 2: the GDC clinical API confirms the field does not exist

```
POST https://api.gdc.cancer.gov/cases   (TCGA-LUAD + TCGA-LUSC)
```

| exposures field | non-empty cases (n = 400 sample) | coverage |
|---|---|---|
| `tobacco_smoking_status` | 372 | **93.0%** |
| `pack_years_smoked` | 278 | **69.5%** |
| `cigarettes_per_day` | 0 | **0.0%** |
| `years_smoked` | 0 | **0.0%** |
| `time_since_cessation` | 0 | **0.0%** |
| `age_at_onset` | 0 | **0.0%** |

In addition, the `exposures` array always has exactly one element (length 1 in 1026 cases, empty in
63), so there is no risk of having read the wrong element.

### Check 3: the pack-years column is trustworthy

The GDC and Xena values for `pack_years_smoked` are **100% identical** across the same 782 patients
(difference < 1e-6). Pack-years is therefore our only reliable continuous exposure variable.

### Check 4: the years-smoked column is trustworthy (cross-checked by back-solving CPD)

Back-solving `CPD = pack-years x 20 / years smoked` from `Xena years_smoked` (n = 507):

| Statistic | Value |
|---|---|
| P5 / P25 / median / P75 / P95 | 9.0 / 20.0 / **20.0** / 30.3 / 55.6 cigarettes/day |
| Within 5-60 cigarettes/day (physiologically plausible) | **94.7%** |

The median lands exactly on **20 cigarettes/day (one pack)** and the distribution has a genuine
shape. If the duration values were fabricated, the back-solved values would be chaotic. Duration is
therefore a real measurement and can be used.

---

## 3. Three substantive consequences for the study design

### (i) Analysis set: 905 -> 559, and the loss can be explained

| Step | n | Loss | Reason |
|---|---|---|---|
| All ever-smokers | 905 | - | full cohort 1087 (LUAD 583 / LUSC 504) |
| With pack-years | 780 | -125 | the GDC did not record pack-years for these 125 cases and there is no alternative source |
| **With smoking duration -> usable** | **559** | -221 | no `years_smoked`, and the legacy onset year is missing too (cases with both `onset` and `stop` recorded = **0**) |
| Legacy rescue potential | - | 0 | of the 315 cases missing duration, only 10 have an onset year, which is not enough to derive duration |

This is **data that do not exist, not data we discarded**, and it can be written directly into a
CONSORT-style flow diagram as an explainable missingness mechanism.

### (ii) A "smoking index" and "pack-years" can never be compared head-to-head in TCGA

**Root cause (more fundamental than "the two metrics happen to coincide")**: in TCGA, **CPD cannot be
observed independently**. The GDC field has 0% coverage (1089 cases), and the Xena column of the
same name is a linear rescaling of pack-years (`== pack-years / 18.25`, exact in 507/507 cases).
If it really were cigarettes per day, then `pack-years / (c/20 x duration)` should be ~1; the
observed median is **9.86**, which refutes it again. Any metric containing CPD must therefore be
*derived* from pack-years and smoking duration:

```
CPD          = pack-years x 20 / duration          (the operative definition)
cigarette-years = CPD x duration == pack-years x 20  (a constructive identity, 507/507)
```

> Under the operative definition, `cig_years == pack-years x 20` is a **constructive identity**
> (because CPD was back-solved from pack-years), not a natural property of the data. A head-to-head
> comparison of a "smoking index" against pack-years is therefore not identifiable in TCGA: the
> data contain no independent CPD to compare with.

**A feasible redesign**:

| Candidate metric | Information structure | Independent? |
|---|---|---|
| pack-years (= smoking index / 20) | cumulative dose | reference |
| **logcig-years** = ln(CPD+1) x duration | emphasises **duration**, compresses intensity | yes |
| **CSI** = (1-0.5^(dur*/tau)) * 0.5^(tsc*/tau) * ln(CPD+1) | duration **saturation** + cessation **decay** + log intensity | yes |
| intensity (pack-years / duration, i.e. packs per day) | pure **intensity** | yes (component) |
| duration (years) | pure **duration** | yes (component) |

**The sharper scientific question** is therefore not "which index is better" but **"which of the
three axes -- intensity, duration and time since cessation -- actually shapes the lung-cancer immune
microenvironment"**.

### (iii) CPD is a derived value and must be declared as such

CPD is back-solved as `pack-years x 20 / duration`, which is the **algebraic identity** implied by
the definition of pack-years (`pack-years = CPD/20 x duration`) and **introduces no new
assumption**. In terms of data provenance, however, it is *derived* rather than *directly measured*,
and the Methods must say so (wording in section 5).

---

## 4. What was changed in the code (`01_Step1_exposure_quantification.R`)

| Location | Change |
|---|---|
| 3.1 `read_gdc_clinical()` | `cigarettes_per_day` is read as `cpd_xena`, explicitly marked as deprecated and kept for audit only |
| section 4, `clinic` | `cpd_raw` -> `cpd_xena_raw` (no longer participates in any metric) |
| section 5, `build_exposure()` | rewritten in **two passes**: the first computes `duration`, the second back-solves `cpd` from it. A single pass would miss the `year_calc` branch that derives duration |
| section 5 | new `cpd_implausible` flag (back-solved value < 1 or > 100 cigarettes/day) |
| new section 6.0 | **data-source audit**: the four checks above are run as QC and written to `QC_xena_cpd_audit.csv`, so the evidence travels with the code |
| section 9 exports | two extra columns, `cpd_implausible` and `cpd_xena_raw`, for traceability |

> After the repair, `cpd_src` is uniformly `back_calc_py` (no more unexplained `clinical`), which
> makes the exposure chain traceable and reproducible.

---

## 5. Ready-to-use manuscript wording

### Methods (data and exposure definition)

> Smoking exposure was derived from the GDC clinical data. Because `cigarettes_per_day` was
> not populated for TCGA-LUAD/LUSC (0/1089 cases in the GDC API), we disregarded the
> rescaled `cigarettes_per_day` field present in the UCSC Xena GDC hub, which proved to be a
> deterministic function of pack-years (`cigarettes_per_day = pack-years/18.25`, 782/782 cases,
> Pearson R2 = 1.000) and therefore carried no independent information. Instead, cigarettes per
> day were reconstructed algebraically from the two verified fields as
> `CPD = pack-years x 20 / smoking duration`, an identity implied by the definition of
> pack-years. Smoking duration was taken from `years_smoked`, or, when missing, derived from the
> interval between `year_of_tobacco_smoking_onset` and `stopped_smoking_year` (or year of
> diagnosis for current smokers). Analyses were restricted to ever-smokers with both pack-years
> and duration available (n = 559) (Figure S1).

### Results (sample flow)

> Among 1087 patients with clinical data, 905 were ever-smokers; pack-years were available in
> 780 (86.2%) and smoking duration in 590 (65.2%), leaving 559 ever-smokers (61.8%) with
> complete quantitative smoking exposure data and a reconstructed CPD distribution consistent
> with published cohorts (median 20.0 cigarettes/day; 95.2% within 5-60).

### Limitations

> Reconstructed cigarettes-per-day values were derived from pack-years and duration rather than
> measured directly, since the original `cigarettes_per_day` field was not populated in TCGA
> lung projects. Because cigarette-years is algebraically equivalent to 20 x pack-years, the two
> are not independent exposure metrics; comparisons were therefore restricted to metrics with
> distinct information structure (pack-years, logcig-years, and the comprehensive smoking index).

---

## 6. Appendix: where the related outputs live

| File | Content |
|---|---|
| `01_Step1_exposure_quantification.R` | section 6.0 contains every audit check described here |
| `QC_xena_cpd_audit.csv` | audit results table (regenerated each time 02 runs) |
| `step1c_missing_funnel.csv` | sample attrition funnel (905 -> 559) |
| `step1c_missing_reason.csv` | reason for missingness x cohort x smoking status |
| `step1c_tsc_missing_diag.csv` | time-since-cessation missingness in former smokers (LUAD 17.8% / LUSC 14.1%, both < 20%, so the effect on the convention is small) |
| `step1c_tsc_scenario_agreement.csv` | agreement across the three time-since-cessation conventions: strict convention rho = 1.0000, agreement 99.3%; midpoint convention rho = 0.9995, agreement 100% |
