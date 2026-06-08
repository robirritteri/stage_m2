#!/usr/bin/env Rscript

# ==========================================================
# Script : pca_bc04.R
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Analyse en composantes principales (ACP/PCA) des profils
#     transcriptomiques de Streptococcus thermophilus (bc04)
#   - Transformation des comptages avec la méthode VST de DESeq2
#   - Visualisation des échantillons selon la condition expérimentale
#     et le stade de fabrication / affinage
#
# Entrées :
#   - Matrice de comptages featureCounts nettoyée pour DESeq2
#     (exemple : bc04.DESeq2.cleaned.tsv)
#   - Répertoire de sortie
#
# Sorties :
#   - PCA_vst.csv
#   - PCA_vst.png
#
# Usage :
#   Rscript pca_bc04.R <counts_tsv> <outdir>
#
# Exemple :
#   Rscript pca_bc04.R \
#     /home/rbirritteri/work/data/counts/bc04.DESeq2.cleaned.tsv \
#     /home/rbirritteri/work/data/pca_bc04
#
# ==========================================================

# ==========================
# Chargement des packages
# ==========================
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ==========================
# Arguments utilisateur
# ==========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage : pca_bc04.R <counts_tsv> <outdir>")
}

counts_tsv <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================
# Fonctions utilitaires
# ==========================

# Ouvre un fichier PNG de manière compatible avec un serveur de calcul
# sans interface graphique.
safe_png <- function(filename, width = 1000, height = 800, res = 150) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPNG(filename, width = width, height = height, res = res)
  } else {
    png(filename, width = width, height = height, res = res, type = "cairo")
  }
}

# Ferme le périphérique graphique.
safe_off <- function() {
  try(dev.off(), silent = TRUE)
}

# Nettoie le nom d'échantillon en supprimant le suffixe de génome bcXX.
clean_sample <- function(x) {
  sub("\\.bc[0-9]+$", "", x)
}

# Extrait le stade expérimental à partir du nom d'échantillon.
stage_from <- function(x) {
  x <- clean_sample(x)
  sub("-[0-9]+$", "", x)
}

# Extrait le numéro de réplicat.
rep_from <- function(x) {
  x <- clean_sample(x)
  as.integer(sub("^.*-([0-9]+)$", "\\1", x))
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

# ==========================
# Lecture des comptages
# ==========================

cts <- read_tsv(counts_tsv, show_col_types = FALSE)

# La première colonne contient les identifiants de gènes.
gene_col <- cts[[1]]

# Les autres colonnes correspondent aux échantillons.
mat <- as.matrix(cts[, -1])
rownames(mat) <- gene_col

samples <- colnames(mat)

# ==========================
# Construction des métadonnées
# ==========================

meta <- data.frame(
  sample = samples,
  stage = stage_from(samples),
  rep = rep_from(samples),
  condition = cond_from_rep(rep_from(samples)),
  row.names = samples,
  stringsAsFactors = FALSE
)

# Définit l'ordre biologique des stades et l'ordre des conditions.
meta$stage <- factor(meta$stage, levels = c("10H", "Fd", "Fp1", "Fp2", "Fp3"))
meta$condition <- factor(meta$condition, levels = c("O", "N"))

# ==========================
# Transformation VST avec DESeq2
# ==========================

# Création de l'objet DESeq2. Le modèle inclut le stade et la condition
# afin que la transformation tienne compte de la structure expérimentale.
dds <- DESeqDataSetFromMatrix(
  countData = round(mat),
  colData = meta,
  design = ~ stage + condition
)

# Préfiltrage des gènes très faiblement exprimés.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Transformation stabilisatrice de variance.
# Elle rend les données plus adaptées à une ACP.
vsd <- vst(dds, blind = FALSE)

# ==========================
# Analyse en composantes principales
# ==========================

# L'ACP est réalisée sur les échantillons à partir des valeurs VST.
pca <- prcomp(t(assay(vsd)))

# Pourcentage de variance expliqué par chaque axe.
pct <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# Table utilisée pour le graphique.
pca_df <- data.frame(pca$x[, 1:2], meta)

# Sauvegarde des coordonnées des échantillons.
write_csv(pca_df, file.path(outdir, "PCA_vst.csv"))

# Représentation graphique de la PCA.
p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = stage)) +
  geom_point(size = 3) +
  labs(
    title = "ACP sur données VST",
    x = sprintf("PC1 (%.1f%%)", pct[1]),
    y = sprintf("PC2 (%.1f%%)", pct[2])
  ) +
  theme_bw()

safe_png(file.path(outdir, "PCA_vst.png"), 1200, 900, 150)
print(p)
safe_off()

message("Analyse terminée : ", outdir)
