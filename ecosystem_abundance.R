#!/usr/bin/env Rscript

# ==========================================================
# Script : ecosystem_abundance.R
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Estimation de l'abondance transcriptionnelle relative des espèces
#     du consortium bactérien à partir des fichiers summary de featureCounts
#   - Utilisation du nombre de lectures assignées comme proxy
#     d'abondance transcriptionnelle
#   - Génération de graphiques par échantillon, par stade et par condition
#
# Entrées :
#   - Répertoire contenant les fichiers *.counts.txt.summary
#   - Répertoire de sortie
#   - Liste des espèces à analyser, séparées par des virgules
#     (exemple : bc01,bc02,bc03,bc04)
#
# Sorties :
#   - relative_abundance_by_sample.csv
#   - relative_abundance_summary_stage_condition.csv
#   - stacked_relative_abundance_by_sample.png
#   - stacked_relative_abundance_mean_stage_condition.png
#
# Usage :
#   Rscript ecosystem_abundance.R <counts_dir> <outdir> <species_csv>
#
# Exemple :
#   Rscript ecosystem_abundance.R \
#     /home/rbirritteri/work/data/counts \
#     /home/rbirritteri/work/data/abundance \
#     bc01,bc02,bc03,bc04
#
# ==========================================================

# ==========================
# Chargement des packages
# ==========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ==========================
# Arguments utilisateur
# ==========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage : ecosystem_abundance.R <counts_dir> <outdir> <species_csv>")
}

counts_dir <- args[1]
outdir <- args[2]
species_list <- strsplit(args[3], ",")[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================
# Fonctions utilitaires
# ==========================

# Ouvre un périphérique PNG robuste, compatible avec les environnements
# de calcul sans interface graphique.
safe_png <- function(filename, width = 1800, height = 700, res = 150) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPNG(filename, width = width, height = height, res = res)
  } else {
    png(filename, width = width, height = height, res = res, type = "cairo")
  }
}

# Ferme proprement le périphérique graphique.
safe_off <- function() {
  try(dev.off(), silent = TRUE)
}

# Nettoie les noms d'échantillons provenant des colonnes featureCounts.
clean_sample <- function(x) {
  x <- sub("^.*/", "", x)
  x <- sub("\\.sorted\\.bam$", "", x)
  x <- sub("\\.bc[0-9][0-9]$", "", x)
  x
}

# Extrait le stade expérimental.
stage_from <- function(s) {
  sub("-\\d+$", "", s)
}

# Extrait le numéro de réplicat.
rep_from <- function(s) {
  suppressWarnings(as.integer(sub("^.*-(\\d+)$", "\\1", s)))
}

# Déduit la condition expérimentale à partir du numéro de réplicat.
cond_from_rep <- function(num) {
  res <- ifelse(
    num <= 4,
    ifelse(num %% 2 == 1, "O", "N"),
    ifelse(num %% 2 == 1, "N", "O")
  )
  res[is.na(num)] <- NA
  return(res)
}

# Lit, pour une espèce donnée, le nombre de lectures assignées dans
# le fichier summary produit par featureCounts.
read_assigned <- function(spec) {
  f <- file.path(counts_dir, paste0(spec, ".counts.txt.summary"))
  if (!file.exists(f)) stop("Fichier manquant : ", f)

  tab <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)

  # La ligne "Assigned" correspond aux lectures assignées à une feature.
  row <- tab %>% filter(Status == "Assigned")
  vals <- as.numeric(row[1, -1])
  cols <- colnames(tab)[-1]
  samples <- vapply(cols, clean_sample, character(1))

  data.frame(sample = samples, assigned = vals, species = spec)
}

# ==========================
# Lecture et mise en forme des données
# ==========================

# Combine les résultats de toutes les espèces demandées.
df <- bind_rows(lapply(species_list, read_assigned)) %>%
  mutate(
    stage = stage_from(sample),
    rep = rep_from(sample),
    condition = cond_from_rep(rep)
  )

# Calcule l'abondance relative par échantillon.
# Pour chaque échantillon, les reads assignés à une espèce sont divisés
# par le total des reads assignés à toutes les espèces.
df_rel <- df %>%
  group_by(sample) %>%
  mutate(
    total_assigned = sum(assigned, na.rm = TRUE),
    rel = ifelse(total_assigned > 0, assigned / total_assigned, NA_real_)
  ) %>%
  ungroup()

df$condition <- factor(df$condition, levels = c("O", "N"))

# Sauvegarde des abondances relatives par échantillon.
write_csv(df_rel, file.path(outdir, "relative_abundance_by_sample.csv"))

# ==========================
# Graphique 1 : abondance relative par échantillon
# ==========================

# Ordonne les échantillons par stade puis par numéro de réplicat.
stage_levels <- c("10H", "Fd", "Fp1", "Fp2", "Fp3")
ord <- df_rel %>%
  distinct(sample, stage, rep) %>%
  arrange(factor(stage, levels = stage_levels), rep, sample) %>%
  pull(sample)

df_rel$sample <- factor(df_rel$sample, levels = ord)
df_rel$stage <- factor(df_rel$stage, levels = stage_levels)

p1 <- ggplot(df_rel, aes(x = sample, y = rel, fill = species)) +
  geom_col(width = 0.9) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Abondance transcriptionnelle relative (proxy = Assigned featureCounts)",
    x = "Échantillon",
    y = "Abondance relative"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing.x = unit(0.2, "lines")
  )

safe_png(file.path(outdir, "stacked_relative_abundance_by_sample.png"), 2400, 700, 150)
print(p1)
safe_off()

# ==========================
# Résumé par stade et condition
# ==========================

df_sum <- df_rel %>%
  group_by(stage, condition, species) %>%
  summarise(
    mean_rel = mean(rel, na.rm = TRUE),
    sd_rel = sd(rel, na.rm = TRUE),
    n = sum(!is.na(rel)),
    .groups = "drop"
  )

write_csv(df_sum, file.path(outdir, "relative_abundance_summary_stage_condition.csv"))

# ==========================
# Graphique 2 : moyenne par stade et condition
# ==========================

p2 <- ggplot(df_sum, aes(x = condition, y = mean_rel, fill = species)) +
  geom_col(width = 0.8) +
  facet_wrap(~stage, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Abondance relative moyenne (stade × condition)",
    x = "Condition (N sans O2 ; O avec O2)",
    y = "Moyenne d'abondance relative"
  ) +
  theme_bw()

safe_png(file.path(outdir, "stacked_relative_abundance_mean_stage_condition.png"), 1700, 500, 150)
print(p2)
safe_off()

message("Analyse terminée : ", outdir)
