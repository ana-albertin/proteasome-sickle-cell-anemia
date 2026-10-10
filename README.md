# Proteasome activity in sickle cell anemia

Proteasome activity (trypsin-like, T-L; chymotrypsin-like, ChT-L), 20S core and
immunoproteasome subunit abundance in healthy controls (CG), untreated patients
(HC-) and patients treated with hydroxycarbamide (hydroxyurea; HC+), all with sickle cell
anemia. HC+ and HC- cells were separated into high-density (HDSC) and low-density (LDSC)
fractions.

## Contents

- `R/`: analysis scripts (Table 1, statistics, figures)
- `data/`: experimental measurements, identified by codes only
- `tables/`: Table 1 (HTML, RTF) and statistical results (`omnibus.csv`, `density.csv`)
- `figures/`: panels A-E (SVG)

Clinical data (`patient_characteristics.csv`: age and blood count) are not included, to
protect patient privacy, and are available on request. `R/01_table1.R` needs that file; its
output is in `tables/`. The other scripts run with the data provided. From the repository
root:

```
Rscript R/01_table1.R   # needs data/patient_characteristics.csv
Rscript R/02_stats.R
Rscript R/03_figures.R
```

Developed with R 4.6.1 (pairwise Wilcoxon p-values with tied values need R >= 4.6.0).
Packages: dplyr, tidyr, gt, lme4, lmerTest, emmeans, car, ggplot2, scales, svglite.

## Data (`data/`)

Patients are identified by codes (`CG01`, `HCpos01`, `HCneg01`, ...). Group codes: `CG`,
`HCpos` (HC+), `HCneg` (HC-). Empty cells are values not measured.

| File | Content |
|---|---|
| `proteasome_activity.csv` | T-L and ChT-L activity (FU/h/µg protein) in total cell extract (`TL_TCE`, `ChTL_TCE`) and after resin purification (`ChTL_resin`). `density` is empty for CG (not fractionated). `excluded = TRUE`: one suspected measurement error, left out of all analyses |
| `wb_20S.csv` | 20S core band density (CG; LDSC of HC+/HC-), normalized to one reference CG sample |
| `wb_immunoproteasome.csv` | β1i and β2i band density, normalized to HC- (n = 1 per group) |
| `patient_characteristics.csv` | Not included: age and blood count per patient (blood count for HC+/HC- only) |

## Analysis

Values are log10-transformed. For each assay:

- **CG vs. HC+ vs. HC-**: one value per patient (HDSC and LDSC averaged; a patient with
  a single fraction contributes that value). The 20S core uses one value per sample (LDSC
  in HC+/HC-). One-way ANOVA with Tukey, or Kruskal-Wallis with pairwise Wilcoxon (BH) if
  the residuals fail Shapiro-Wilk (p < 0.05) (`tables/omnibus.csv`).
- **HDSC vs. LDSC in HC+/HC-**: linear mixed model with patient as random intercept,
  since both fractions come from the same person (`tables/density.csv`, log10 scale).

Pairwise p values are adjusted within each assay (Tukey or BH).
Immunoproteasome subunits are descriptive only.
