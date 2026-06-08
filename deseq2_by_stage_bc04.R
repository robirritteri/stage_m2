#!/usr/bin/env Rscript

# ==========================================================
# Script : deseq2_by_stage_bc04.R
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Analyse différentielle d'expression par stade expérimental
#     pour Streptococcus thermophilus (bc04)
#   - Comparaison de la condition N par rapport à la condition O
#     à chaque stade : 10H, Fd, Fp1, Fp2 et Fp3
#   - Génération des matrices de comptages normalisés, des résultats
#     DESeq2 et d'un tableau récapitulatif par stade
#
# Entrées :
#   - Matrice de comptages featureCounts nettoyée pour DESeq2
#     (exemple : bc04.DESeq2.cleaned.tsv)
#   - Répertoire de sortie
#
# Sorties :
#   - samplesheet_all.tsv
#   - samplesheet_<stade>.tsv
#   - normalized_counts_<stade>.tsv
#   - DE_<stade>_condition_N_vs_O.csv
#   - DE_summary_by_stage.csv
#
# Usage :
#   Rscript deseq2_by_stage_bc04.R <counts_cleaned.tsv> <outdir>
#
# Exemple :
#   Rscript deseq2_by_stage_bc04.R \
#     /home/rbirritteri/work/data/counts/bc04.DESeq2.cleaned.tsv \
#     /home/rbirritteri/work/data/deseq2/deseq2_results_bc04_by_stage
#
# ==========================================================

# Fixe la locale en anglais/C pour limiter les problèmes de format
# numérique ou de messages selon l'environnement d'exécution.
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)

# ==========================
# Chargement des packages
# ==========================
suppressPackageStartupMessages({
  library(DESeq2)
})

# ==========================
# Arguments utilisateur
# ==========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage : deseq2_by_stage_bc04.R <counts_cleaned.tsv> <outdir>\n",
    "Exemple : deseq2_by_stage_bc04.R /home/rbirritteri/work/data/counts/bc04.DESeq2.cleaned.tsv /home/rbirritteri/work/data/deseq2/deseq2_results_bc04_by_stage\n"
  )
}

counts_file <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Fichier de comptages : ", counts_file)
message("Répertoire de sortie : ", out_dir)

# ==========================
# Fonctions utilitaires
# ==========================

# Nettoie les noms de colonnes issus de featureCounts / BAM afin de récupérer
# uniquement le nom d'échantillon.
parse_sample <- function(x) {
  x0 <- x
  x0 <- sub("\\.sorted.*$", "", x0)
  x0 <- sub("\\.bam$", "", x0)
  x0 <- sub("\\.bc\\d\\d(\\..*)?$", "", x0)
  x0 <- sub("\\.counts.*$", "", x0)
  x0
}

# Extrait le stade expérimental à partir du nom d'échantillon.
# Exemple : Fp1-03 -> Fp1.
stage_from <- function(sample) {
  stg <- sub("-\\d+$", "", sample)
  stg[stg == sample] <- NA_character_
  stg
}

# Extrait le numéro de réplicat à partir du nom d'échantillon.
# Exemple : Fp1-03 -> 3.
rep_from <- function(sample) {
  n <- suppressWarnings(as.integer(sub("^.*-(\\d+)$", "\\1", sample)))
  ifelse(is.na(n) | is.nan(n), NA_integer_, n)
}

# Déduit la condition expérimentale à partir du numéro de réplicat.
# Les réplicats 1 à 4 et 5 à 8 suivent une alternance O/N inversée.
cond_from_rep <- function(num) {
  res <- ifelse(
    num <= 4,
    ifelse(num %% 2 == 1, "O", "N"),
    ifelse(num %% 2 == 1, "N", "O")
  )
  res[is.na(num)] <- NA_character_
  res
}

# Ordre attendu des stades expérimentaux.
stage_levels <- c("10H", "Fd", "Fp1", "Fp2", "Fp3")

# ==========================
# Lecture de la matrice de comptages
# ==========================
tab <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)

# La première colonne correspond aux identifiants de gènes.
gene_id <- tab[[1]]

# Les autres colonnes correspondent aux échantillons.
count_mat <- as.matrix(tab[, -1, drop = FALSE])
rownames(count_mat) <- gene_id

# ==========================
# Construction des métadonnées
# ==========================
samples_raw <- colnames(count_mat)
samples <- vapply(samples_raw, parse_sample, character(1))

stage <- stage_from(samples)
repn  <- rep_from(samples)
condition <- cond_from_rep(repn)

# Vérifie que le stade et la condition ont bien été identifiés pour
# chaque échantillon.
bad <- which(is.na(stage) | is.na(condition))
if (length(bad) > 0) {
  msg <- paste0(
    "Erreur de métadonnées : stade ou condition non identifié pour ces colonnes :\n",
    paste0(
      " - colonne='", samples_raw[bad],
      "' échantillon='", samples[bad],
      "' stade='", stage[bad],
      "' rep='", repn[bad],
      "' condition='", condition[bad],
      "'", collapse = "\n"
    )
  )
  stop(msg)
}

# Vérifie que les stades détectés correspondent aux stades attendus.
unknown_stage <- unique(stage[!(stage %in% stage_levels)])
if (length(unknown_stage) > 0) {
  stop(
    "Stades inconnus détectés : ", paste(unknown_stage, collapse = ", "),
    "\nStades attendus : ", paste(stage_levels, collapse = ", ")
  )
}

# Table de métadonnées utilisée par DESeq2.
coldata <- data.frame(
  sample = samples,
  stage = factor(stage, levels = stage_levels),
  condition = factor(condition, levels = c("O", "N")),
  rep = repn,
  stringsAsFactors = FALSE
)
rownames(coldata) <- samples_raw

# Sauvegarde de la table complète des métadonnées.
write.table(
  coldata,
  file = file.path(out_dir, "samplesheet_all.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)

# ==========================
# Préfiltrage global des gènes
# ==========================
# Les gènes avec très peu de lectures sur l'ensemble des échantillons sont
# retirés afin de limiter le bruit statistique.
keep_genes <- rowSums(count_mat) >= 10
count_mat_f <- count_mat[keep_genes, , drop = FALSE]
message("Gènes conservés après préfiltrage global : ", nrow(count_mat_f), " / ", nrow(count_mat))

# ==========================
# Analyse DESeq2 stade par stade
# ==========================
summary_list <- list()

for (stg in stage_levels) {
  message("===== Stade : ", stg, " =====")

  stage_dir <- file.path(out_dir, stg)
  dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)

  # Sélection des échantillons correspondant au stade courant.
  idx <- which(coldata$stage == stg)
  if (length(idx) == 0) {
    message("Aucun échantillon pour le stade ", stg, " : étape ignorée.")
    next
  }

  coldata_stg <- droplevels(coldata[idx, , drop = FALSE])
  count_mat_stg <- count_mat_f[, rownames(coldata_stg), drop = FALSE]

  # Vérifie le nombre d'échantillons par condition.
  cond_tab <- table(coldata_stg$condition)
  message("Nombre d'échantillons par condition pour ", stg, " : ",
          paste(names(cond_tab), cond_tab, collapse = ", "))

  if (!all(c("O", "N") %in% names(cond_tab))) {
    message("Une condition est absente au stade ", stg, " : étape ignorée.")
    next
  }

  if (any(cond_tab < 2)) {
    message("Attention : moins de 2 réplicats dans une condition pour le stade ", stg)
  }

  # Préfiltrage spécifique au stade analysé.
  keep_stage <- rowSums(count_mat_stg) >= 10
  count_mat_stg_f <- count_mat_stg[keep_stage, , drop = FALSE]
  message("Gènes conservés pour ", stg, " : ", nrow(count_mat_stg_f), " / ", nrow(count_mat_stg))

  # Sauvegarde des métadonnées du stade.
  write.table(
    coldata_stg,
    file = file.path(stage_dir, paste0("samplesheet_", stg, ".tsv")),
    sep = "\t", quote = FALSE, row.names = TRUE
  )

  # Création de l'objet DESeq2.
  # Le modèle teste l'effet de la condition au sein d'un stade donné.
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat_stg_f),
    colData = coldata_stg,
    design = ~ condition
  )

  # Normalisation et estimation du modèle statistique.
  dds <- DESeq(dds)

  # Extraction et sauvegarde des comptages normalisés.
  norm_counts <- counts(dds, normalized = TRUE)
  write.table(
    norm_counts,
    file = file.path(stage_dir, paste0("normalized_counts_", stg, ".tsv")),
    sep = "\t", quote = FALSE
  )

  # Résultats de l'analyse différentielle.
  # Le contraste N vs O indique les gènes plus ou moins exprimés en condition N
  # par rapport à la condition O.
  res <- results(dds, contrast = c("condition", "N", "O"))
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj), ]

  out_csv <- file.path(stage_dir, paste0("DE_", stg, "_condition_N_vs_O.csv"))
  write.csv(res_df, file = out_csv, row.names = TRUE)

  # Tableau récapitulatif du nombre de gènes testés et significatifs.
  summary_list[[stg]] <- data.frame(
    stage = stg,
    n_samples = nrow(coldata_stg),
    n_O = unname(cond_tab["O"]),
    n_N = unname(cond_tab["N"]),
    n_genes_tested = nrow(res_df),
    n_padj_0.05 = sum(res_df$padj < 0.05, na.rm = TRUE),
    n_padj_0.10 = sum(res_df$padj < 0.10, na.rm = TRUE),
    n_padj_0.05_lfc0.5 = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0.5, na.rm = TRUE),
    n_padj_0.05_lfc1 = sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, na.rm = TRUE)
  )

  message("Résultats écrits : ", out_csv)
}

# ==========================
# Résumé final
# ==========================
if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.csv(
    summary_df,
    file = file.path(out_dir, "DE_summary_by_stage.csv"),
    row.names = FALSE
  )
  message("Résumé écrit : ", file.path(out_dir, "DE_summary_by_stage.csv"))
}

message("Analyse terminée.")
