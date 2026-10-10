# Table 1: patient characteristics (HC+ vs. HC-), exported as HTML and RTF.

library(dplyr)
library(tidyr)
library(gt)

dir.create("tables", showWarnings = FALSE)

pc <- read.csv("data/patient_characteristics.csv")

# Reference ranges are standard adult values, not computed from the data.
# dec = decimals shown.
vars <- tribble(
  ~variable,                   ~label,         ~unit,                 ~dec, ~ref,
  "age",                       "Age",          "years",               1,    "-",
  "erythrocytes_millions_mm3", "Erythrocytes", "millions/mm\u00b3",   2,    "4.0-5.0",
  "hematocrit_pct",            "Hematocrit",   "%",                   1,    "36.0-47.0",
  "hemoglobin_g_dl",           "Hemoglobin",   "g/dL",                1,    "12.0-16.0",
  "mcv_fl",                    "MCV",          "fL",                  1,    "80.0-100.0",
  "mch_pg",                    "MCH",          "pg",                  1,    "27.0-31.0",
  "mchc_g_dl",                 "MCHC",         "g/dL",                1,    "32.0-36.0",
  "rdw_pct",                   "RDW",          "%",                   1,    "12.0-16.0",
  "leukocytes_per_mm3",        "Leukocytes",   "cells/mm\u00b3",      0,    "4,000-11,000",
  "granulocytes_per_mm3",      "Granulocytes", "cells/mm\u00b3",      0,    "-",
  "lymphocytes_per_mm3",       "Lymphocytes",  "cells/mm\u00b3",      0,    "-",
  "monocytes_per_mm3",         "Monocytes",    "cells/mm\u00b3",      0,    "-",
  "platelets_k_ul",            "Platelets",    "x10\u00b3/\u00b5L",  0,    "-"
)

summ <- pc %>%
  pivot_longer(-c(patient_id, group), names_to = "variable") %>%
  filter(!is.na(value)) %>%
  left_join(select(vars, variable, dec), by = "variable") %>%
  group_by(group, variable) %>%
  summarize(cell = sprintf("%s \u00b1 %s (n=%d)",
                           formatC(mean(value), format = "f", digits = dec[1]),
                           formatC(sd(value), format = "f", digits = dec[1]), n()),
            .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = cell)

tbl <- vars %>%
  left_join(summ, by = "variable") %>%
  select(label, unit, CG, HCpos, HCneg, ref) %>%
  gt(id = "table1") %>%
  tab_header(title = "Table 1. Patient characteristics") %>%
  cols_label(label = "Variable", unit = "Unit", CG = "CG (mean \u00b1 SD)",
             HCpos = "HC+ (mean \u00b1 SD)", HCneg = "HC- (mean \u00b1 SD)",
             ref = "Reference range") %>%
  cols_align("center", columns = c(unit, CG, HCpos, HCneg, ref)) %>%
  sub_missing(missing_text = "-") %>%
  tab_source_note("CG: healthy controls (blood count not collected); HC+: treated with hydroxycarbamide (hydroxyurea); HC-: untreated.") %>%
  tab_style(cell_text(weight = "bold"), cells_title()) %>%
  tab_style(cell_text(weight = "bold"), cells_column_labels()) %>%
  opt_table_font("Arial")

gtsave(tbl, "tables/table1.html")
gtsave(tbl, "tables/table1.rtf")
