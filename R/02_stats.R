# Statistics for proteasome activity (T-L, ChT-L) and 20S core abundance.
# All models use log10 values.
#
# omnibus: CG vs. HC+ vs. HC-, one value per patient (HDSC/LDSC averaged for
#   HC+/HC-; patients with a single fraction contribute that value). One-way
#   ANOVA + Tukey if the residuals pass Shapiro-Wilk (p >= 0.05); otherwise
#   Kruskal-Wallis + pairwise Wilcoxon (BH). Levene's p is reported only.
#   Pairwise p-values are adjusted within each assay (Tukey/BH).
# density: HC+/HC- only. Linear mixed model with patient as random intercept
#   (HDSC and LDSC come from the same person); HDSC vs. LDSC within each group.

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

dir.create("tables", showWarnings = FALSE)

omnibus <- function(d) {  # d: y (log10 value), group
  d$group <- factor(d$group, levels = c("CG", "HCpos", "HCneg"))
  fit <- lm(y ~ group, data = d)
  shapiro <- shapiro.test(residuals(fit))$p.value
  if (shapiro >= 0.05) {
    test <- "One-way ANOVA"
    p <- anova(fit)[["Pr(>F)"]][1]
    pw <- as.data.frame(pairs(emmeans(fit, ~ group)))$p.value  # Tukey
  } else {
    test <- "Kruskal-Wallis"
    p <- kruskal.test(y ~ group, data = d)$p.value
    m <- pairwise.wilcox.test(d$y, d$group, p.adjust.method = "BH")$p.value
    pw <- c(m["HCpos", "CG"], m["HCneg", "CG"], m["HCneg", "HCpos"])
  }
  data.frame(test = test, shapiro_p = shapiro,
             levene_p = car::leveneTest(y ~ group, data = d)[1, "Pr(>F)"],
             omnibus_p = p, p_cg_hcpos = pw[1], p_cg_hcneg = pw[2], p_hcpos_hcneg = pw[3])
}

density_test <- function(d) {  # d: HC+/HC- rows with y, group, density, patient_id
  d$group <- factor(d$group, levels = c("HCpos", "HCneg"))
  fit <- lmer(y ~ group * density + (1 | patient_id), data = d)
  as.data.frame(emmeans(fit, pairwise ~ density | group, lmer.df = "satterthwaite")$contrasts) %>%
    mutate(interaction_p = anova(fit)["group:density", "Pr(>F)"],
           shapiro_p = shapiro.test(residuals(fit))$p.value)
}

act <- read.csv("data/proteasome_activity.csv", na.strings = "") %>%
  filter(!excluded) %>%
  mutate(y = log10(activity))
wb <- read.csv("data/wb_20S.csv") %>% mutate(y = log10(normalized))
assays <- c("TL_TCE", "ChTL_TCE", "ChTL_resin")

omni <- bind_rows(
  lapply(assays, function(a) {
    act %>% filter(assay == a) %>%
      group_by(patient_id, group) %>% summarize(y = mean(y), .groups = "drop") %>%
      omnibus() %>% mutate(assay = a)
  }),
  omnibus(wb) %>% mutate(assay = "20S")
) %>% relocate(assay)

dens <- bind_rows(lapply(assays, function(a) {
  density_test(filter(act, assay == a, group != "CG")) %>% mutate(assay = a)
})) %>% relocate(assay)

write.csv(omni, "tables/omnibus.csv", row.names = FALSE)
write.csv(dens, "tables/density.csv", row.names = FALSE)
print(omni)
print(dens)
