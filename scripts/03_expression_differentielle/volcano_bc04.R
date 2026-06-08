#!/usr/bin/env Rscript

# ==========================================================
# Script : volcano_bc04.R
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Génération de volcano plots à partir des résultats DESeq2
#     obtenus pour Streptococcus thermophilus (bc04)
#   - Un graphique est produit pour chaque stade expérimental
#   - Les gènes sont classés comme significatifs selon un seuil
#     de p-value ajustée et un seuil de log2FoldChange
#
# Entrées :
#   - Répertoire contenant les fichiers DE_<stade>_N_vs_O.csv
#     ou DE_<stade>_condition_N_vs_O.csv
#   - Répertoire de sortie
#   - Seuil optionnel de p-value ajustée
#   - Seuil optionnel de log2FoldChange
#
# Sorties :
#   - volcano_bc04_<stade>.png
#
# Usage :
#   Rscript volcano_bc04.R <deseq_dir> <outdir> [padj] [lfc]
#
# Exemple :
#   Rscript volcano_bc04.R \
#     /home/rbirritteri/work/data/deseq2/deseq2_results_bc04_by_stage \
#     /home/rbirritteri/work/data/volcano_bc04 \
#     0.05 0.5
#
# ==========================================================

# ==========================
# Chargement des packages
# ==========================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ==========================
# Arguments utilisateur
# ==========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage : volcano_bc04.R <deseq_dir> <outdir> [padj] [lfc]")
}

deseq_dir <- args[1]
outdir <- args[2]

# Seuils utilisés pour définir les gènes significatifs.
# Par défaut : padj < 0,05 et |log2FoldChange| > 0,5.
padj_thr <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)
lfc_thr  <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.5)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================
# Fonctions utilitaires
# ==========================

# Ouvre un fichier PNG compatible avec un environnement de calcul sans écran.
safe_png <- function(filename, width = 900, height = 800, res = 150) {
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

# ==========================
# Recherche des fichiers DESeq2
# ==========================

# Les fichiers attendus correspondent aux résultats de comparaison N vs O.
# Deux formats de nom sont acceptés pour rester compatible avec les sorties
# possibles du script DESeq2.
files <- list.files(
  deseq_dir,
  pattern = "^DE_.*(_condition)?_N_vs_O\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(files) == 0) {
  stop("Aucun fichier DE_*_N_vs_O.csv trouvé dans : ", deseq_dir)
}

# Extrait le stade expérimental à partir du nom du fichier.
get_stage <- function(f) {
  x <- basename(f)
  x <- sub("^DE_", "", x)
  x <- sub("_condition_N_vs_O\\.csv$", "", x)
  x <- sub("_N_vs_O\\.csv$", "", x)
  x
}

# ==========================
# Génération des volcano plots
# ==========================

for (f in files) {
  stg <- get_stage(f)
  df <- read_csv(f, show_col_types = FALSE)

  # La première colonne correspond à l'identifiant du gène.
  names(df)[1] <- "gene"

  # Ajoute une colonne indiquant si le gène passe les seuils choisis.
  # neglog10 correspond à -log10(padj), utilisé sur l'axe Y.
  df <- df %>%
    mutate(
      sig = !is.na(padj) & padj < padj_thr & abs(log2FoldChange) > lfc_thr,
      neglog10 = ifelse(!is.na(padj) & padj > 0, -log10(padj), NA_real_)
    )

  # Volcano plot : effet d'expression en X, significativité en Y.
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10, color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    labs(
      title = paste0("Volcano plot bc04 - ", stg, " (N vs O)"),
      x = "log2FoldChange",
      y = "-log10(padj)",
      color = "Significatif"
    ) +
    theme_bw()

  out <- file.path(outdir, paste0("volcano_bc04_", stg, ".png"))
  safe_png(out, 900, 800, 150)
  print(p)
  safe_off()
}

message("Analyse terminée : ", outdir)
