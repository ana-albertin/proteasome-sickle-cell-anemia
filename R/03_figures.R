# Figure panels (SVG): A-C proteasome activity, D 20S core, E immunoproteasome.
# Boxplots (outline only) with the geometric mean marked as "+"; points are
# single patients, joined across HDSC/LDSC. Statistics come from tables/
# (02_stats.R).

library(dplyr)
library(ggplot2)
library(scales)

dir.create("figures", showWarnings = FALSE)

pal <- c(CG = "#A5CCE9", `HC+` = "#A5E9D8", `HC-` = "#E9C1A5")  # colorblind-safe
shapes <- c(HDSC = 24, LDSC = 21, Total = 23)

as_group <- function(g) factor(recode(g, HCpos = "HC+", HCneg = "HC-"), levels = names(pal))
stars <- function(p) as.character(cut(p, c(-Inf, .001, .01, .05, Inf), c("***", "**", "*", "ns"), right = FALSE))
fmt_p <- function(p) ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))
star_size <- 14 / .pt  # significance marks: 14 pt

theme_paper <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.title.position = "plot",
    plot.subtitle = element_text(size = 12, color = "gray25", lineheight = 1.15),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(10, 16, 6, 6)
  )

act  <- read.csv("data/proteasome_activity.csv", na.strings = "") %>% filter(!excluded)
wb   <- read.csv("data/wb_20S.csv")
subu <- read.csv("data/wb_immunoproteasome.csv")
omni <- read.csv("tables/omnibus.csv")
dens <- read.csv("tables/density.csv")

# A-C: activity by group and density; bracket = HDSC vs. LDSC (mixed model).
plot_activity <- function(a, title, ylab) {
  d <- act %>% filter(assay == a) %>%
    mutate(group_lb = as_group(group),
           density = factor(coalesce(density, "Total"), levels = names(shapes)))
  o <- filter(omni, assay == a)
  gm <- d %>% group_by(group_lb, density) %>%
    summarize(m = 10^mean(log10(activity)), .groups = "drop")  # geometric mean
  br <- d %>% filter(group != "CG") %>% group_by(group_lb) %>%
    summarize(y = max(activity) * 1.3, .groups = "drop") %>%
    left_join(dens %>% filter(assay == a) %>% transmute(group_lb = as_group(group), label = stars(p.value)),
              by = "group_lb")
  subtitle <- sprintf("%s (CG vs. HC+ vs. HC-): %s %s\nCG vs. HC+ %s, CG vs. HC- %s, HC+ vs. HC- %s",
                      o$test, fmt_p(o$omnibus_p), stars(o$omnibus_p),
                      stars(o$p_cg_hcpos), stars(o$p_cg_hcneg), stars(o$p_hcpos_hcneg))

  ggplot(d, aes(density, activity)) +
    geom_boxplot(aes(group = density), fill = NA, color = "gray35", width = 0.5,
                 outlier.shape = NA, linewidth = 0.4) +
    geom_line(data = filter(d, group != "CG"), aes(group = patient_id), color = "gray70", linewidth = 0.3) +
    geom_jitter(aes(fill = group_lb, shape = density), color = "black",
                position = position_jitter(width = 0.08, height = 0, seed = 1), size = 2.3, stroke = 0.3, alpha = 0.85) +
    geom_point(data = gm, aes(y = m), shape = 3, size = 3.2, stroke = 1.3) +
    geom_segment(data = br, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = br, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, vjust = -0.3, size = star_size) +
    facet_grid(. ~ group_lb, scales = "free_x", space = "free_x") +
    scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0.05, 0.16))) +
    scale_fill_manual(values = pal, name = "Group") +
    scale_shape_manual(values = shapes, name = "Density") +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
           shape = guide_legend(order = 2)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_paper
}

# D: 20S core; brackets = pairwise comparisons from the omnibus test.
plot_20S <- function() {
  d <- wb %>% mutate(group_lb = as_group(group))
  o <- filter(omni, assay == "20S")
  gm <- d %>% group_by(group_lb) %>% summarize(m = 10^mean(log10(normalized)), .groups = "drop")
  br <- data.frame(x1 = c(1, 2, 1), x2 = c(2, 3, 3),
                   y = max(d$normalized) * 10^c(0.05, 0.12, 0.24),
                   label = stars(c(o$p_cg_hcpos, o$p_hcpos_hcneg, o$p_cg_hcneg)))

  ggplot(d, aes(group_lb, normalized)) +
    geom_boxplot(aes(group = group_lb), fill = NA, color = "gray35", width = 0.5,
                 outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(aes(fill = group_lb), shape = 21, color = "black",
                position = position_jitter(width = 0.08, height = 0, seed = 1), size = 2.3, stroke = 0.3, alpha = 0.85) +
    geom_point(data = gm, aes(y = m), shape = 3, size = 3.2, stroke = 1.3) +
    geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = br, aes(x = (x1 + x2) / 2, y = y, label = label), inherit.aes = FALSE, vjust = -0.3, size = star_size) +
    scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0.05, 0.12))) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = "20S core",
         subtitle = sprintf("%s (CG vs. HC+ vs. HC-): %s %s", o$test, fmt_p(o$omnibus_p), stars(o$omnibus_p)),
         x = NULL, y = "Optical Densitometry (RU)") +
    theme_paper
}

# E: immunoproteasome subunits (n = 1 per group, descriptive only).
plot_subunits <- function() {
  subu %>%
    mutate(group_lb = as_group(group), subunit = recode(subunit, beta1i = "\u03b21i", beta2i = "\u03b22i")) %>%
    ggplot(aes(subunit, normalized, fill = group_lb)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = pal, name = "Group") +
    labs(title = "Immunoproteasome subunits", subtitle = "n=1/group, descriptive only",
         x = NULL, y = "Relative abundance\n(vs. HC-)") +
    theme_paper +
    theme(legend.position = "right")  # saves height
}

unit <- "(FU/h/\u00b5g protein)"
panels <- list(
  A = plot_activity("ChTL_TCE",   "ChT-L activity (TCE)",   paste("ChT-L activity", unit)),
  B = plot_activity("TL_TCE",     "T-L activity (TCE)",     paste("T-L activity", unit)),
  C = plot_activity("ChTL_resin", "ChT-L activity (resin)", paste("ChT-L activity", unit)),
  D = plot_20S(),
  E = plot_subunits()
)

for (n in names(panels)) {
  h <- if (n == "E") 4.6 / 2 else 4.6  # panel E: half height
  ggsave(sprintf("figures/panel_%s.svg", n), panels[[n]], width = 4.8, height = h, device = "svg")
}
