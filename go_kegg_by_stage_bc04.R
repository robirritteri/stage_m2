#!/usr/bin/env Rscript

# ==========================================================
# Script : go_kegg_by_stage_bc04.R
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Analyse d'enrichissement fonctionnel KEGG par stade expérimental
#     pour Streptococcus thermophilus (bc04)
#   - Utilisation des résultats DESeq2 obtenus pour la comparaison N vs O
#   - Croisement des gènes différentiellement exprimés avec les annotations
#     fonctionnelles eggNOG / KEGG
#   - Réalisation d'une analyse d'enrichissement ORA avec clusterProfiler
#     séparément pour les gènes surexprimés et sous-exprimés en condition N
#
# Entrées :
#   - Répertoire contenant les résultats DESeq2 par stade
#     (DE_<stade>_condition_N_vs_O.csv)
#   - Tableau d'annotation fonctionnelle issu de eggNOG
#     (bc04_annotation.tsv)
#   - Répertoire de sortie
#
# Sorties :
#   - Liste des pathways KEGG utilisés
#   - Tables TERM2GENE et TERM2NAME KEGG
#   - Fichiers de diagnostic par stade
#   - Gènes différentiellement exprimés annotés
#   - Résultats d'enrichissement KEGG pour les gènes up/down
#   - Graphiques dotplot et barplot des pathways enrichis
#   - Tableau résumé des enrichissements par stade
#
# Usage :
#   Rscript go_kegg_by_stage_bc04.R <deseq_stage_dir> <annotation_tsv> <outdir>
#
# Exemple :
#   Rscript go_kegg_by_stage_bc04.R \
#     /home/rbirritteri/work/data/deseq2/deseq2_results_bc04_by_stage \
#     /home/rbirritteri/work/data/functional_annotation/bc04_annotation.tsv \
#     /home/rbirritteri/work/data/kegg/bc04_by_stage
#
# ==========================================================
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)

# ==========================
# Chargement des packages
# ==========================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

# ==========================
# Arguments utilisateur
# ==========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage : Rscript go_kegg_by_stage_bc04.R <deseq_stage_dir> <annotation_tsv> <outdir>\n",
    "Exemple :\n",
    "Rscript go_kegg_by_stage_bc04.R ",
    "/home/rbirritteri/work/data/deseq2/deseq2_results_bc04_by_stage ",
    "/home/rbirritteri/work/data/functional_annotation/bc04_annotation.tsv ",
    "/home/rbirritteri/work/data/kegg/bc04_by_stage\n"
  )
}

deseq_dir <- args[1]
annot_file <- args[2]
out_dir <- args[3]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Ordre biologique des stades de fabrication / affinage.
stage_levels <- c("10H", "Fd", "Fp1", "Fp2", "Fp3")

message("Répertoire DESeq2 : ", deseq_dir)
message("Annotation        : ", annot_file)
message("Répertoire sortie : ", out_dir)

# ==========================
# 1. Lecture des annotations eggNOG
# ==========================

# Le tableau d'annotation est obtenu après parsing du fichier eggNOG-mapper.
# Une ligne correspond à un locus_tag et contient notamment les voies KEGG.
ann <- read_tsv(annot_file, show_col_types = FALSE) %>%
  distinct(locus_tag, .keep_all = TRUE)

# Les colonnes locus_tag et KEGG_Pathway sont indispensables pour relier
# les gènes DESeq2 aux voies fonctionnelles KEGG.
required_cols <- c("locus_tag", "KEGG_Pathway")
missing_cols <- setdiff(required_cols, colnames(ann))

if (length(missing_cols) > 0) {
  stop("Colonnes manquantes : ", paste(missing_cols, collapse = ", "))
}

# ==========================
# 2. Récupération des noms officiels KEGG
# ==========================

# Récupère la correspondance officielle entre identifiants KEGG pathway
# et noms de pathways. Cette table est utilisée pour rendre les sorties
# d'enrichissement plus lisibles.
get_kegg_term2name <- function() {
  url <- "https://rest.kegg.jp/list/pathway/ko"

  tryCatch({
    df <- read.delim(url, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c("KEGG_Pathway", "Description")

    df %>%
      mutate(
        KEGG_Pathway = str_remove(KEGG_Pathway, "^path:"),
        Description = str_remove(Description, " - Reference pathway$")
      ) %>%
      filter(str_detect(KEGG_Pathway, "^ko\\d{5}$")) %>%
      distinct(KEGG_Pathway, .keep_all = TRUE)
  }, error = function(e) {
    message("Impossible de récupérer KEGG TERM2NAME : ", e$message)
    return(NULL)
  })
}

kegg_term2name_all <- get_kegg_term2name()

# ==========================
# 3. Sélection des pathways KEGG de type Metabolism
# ==========================

# Récupère automatiquement les identifiants des pathways classés dans la
# catégorie Metabolism de KEGG BRITE. Cela évite de fournir une liste manuelle.
get_kegg_metabolism_ids <- function() {
  url <- "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext"

  lines <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) {
      stop("Impossible de récupérer KEGG BRITE br08901 : ", e$message)
    }
  )

  in_metabolism <- FALSE
  ids <- character()

  for (line in lines) {

    if (str_detect(line, "^A<b>Metabolism</b>")) {
      in_metabolism <- TRUE
      next
    }

    if (str_detect(line, "^A<b>") && !str_detect(line, "^A<b>Metabolism</b>")) {
      in_metabolism <- FALSE
    }

    if (in_metabolism) {
      hit <- str_match(line, "\\b([0-9]{5})\\b")
      if (!is.na(hit[, 2])) {
        ids <- c(ids, paste0("ko", hit[, 2]))
      }
    }
  }

  unique(ids)
}

metabolism_ids <- get_kegg_metabolism_ids()

if (length(metabolism_ids) == 0) {
  stop("Aucun pathway KEGG Metabolism récupéré.")
}

# Les pathways très globaux sont retirés pour éviter des enrichissements
# trop généraux et peu informatifs biologiquement.
global_metabolism_ids <- c(
  "ko01100", # Metabolic pathways
  "ko01110", # Biosynthesis of secondary metabolites
  "ko01120", # Microbial metabolism in diverse environments
  "ko01200", # Carbon metabolism
  "ko01210", # 2-Oxocarboxylic acid metabolism
  "ko01230"  # Biosynthesis of amino acids
)

metabolism_ids_filtered <- setdiff(metabolism_ids, global_metabolism_ids)

write_csv(
  tibble(KEGG_Pathway = metabolism_ids_filtered),
  file.path(out_dir, "KEGG_metabolism_pathways_used.csv")
)

# ==========================
# 4. Construction des tables TERM2GENE / TERM2NAME
# ==========================

# clusterProfiler utilise une table TERM2GENE associant chaque pathway KEGG
# aux gènes annotés dans ce pathway. Cette table est construite à partir
# des annotations eggNOG de bc04.
ann_kegg_raw <- ann %>%
  select(locus_tag, KEGG_Pathway) %>%
  filter(
    !is.na(KEGG_Pathway),
    KEGG_Pathway != "",
    KEGG_Pathway != "-",
    KEGG_Pathway != "--"
  ) %>%
  mutate(KEGG_Pathway = str_replace_all(KEGG_Pathway, ";", ",")) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(KEGG_Pathway = str_trim(KEGG_Pathway)) %>%
  filter(str_detect(KEGG_Pathway, "^ko\\d{5}$")) %>%
  distinct(KEGG_Pathway, locus_tag)

# Première sélection limitée aux pathways métaboliques retenus.
ann_kegg <- ann_kegg_raw %>%
  filter(KEGG_Pathway %in% metabolism_ids_filtered)

if (!is.null(kegg_term2name_all)) {
  kegg_term2name <- kegg_term2name_all %>%
    filter(KEGG_Pathway %in% unique(ann_kegg$KEGG_Pathway))
} else {
  kegg_term2name <- NULL
}

# Quelques pathways non strictement métaboliques sont conservés car ils peuvent
# être informatifs dans un contexte transcriptomique bactérien.
extra_kegg <- c(
  "ko02010", # ABC transporters
  "ko02020", # Two-component system
  "ko02024", # Quorum sensing
  "ko03010", # Ribosome
  "ko03020", # RNA polymerase
  "ko03030", # DNA replication
  "ko03410", # Base excision repair
  "ko03420"  # Nucleotide excision repair
)

# Table finale utilisée pour l'enrichissement.
ann_kegg <- ann_kegg_raw %>%
  filter(
    KEGG_Pathway %in% metabolism_ids_filtered |
    KEGG_Pathway %in% extra_kegg
  )

write_csv(ann_kegg_raw, file.path(out_dir, "KEGG_TERM2GENE_raw.csv"))
write_csv(ann_kegg, file.path(out_dir, "KEGG_TERM2GENE_metabolism_only.csv"))

if (!is.null(kegg_term2name)) {
  write_csv(kegg_term2name, file.path(out_dir, "KEGG_TERM2NAME_metabolism_only.csv"))
}

message("Pathways KEGG bruts           : ", length(unique(ann_kegg_raw$KEGG_Pathway)))
message("Pathways KEGG metabolism gardés: ", length(unique(ann_kegg$KEGG_Pathway)))
message("Gènes avec KEGG metabolism    : ", length(unique(ann_kegg$locus_tag)))

# ==========================
# 5. Fonction d'enrichissement KEGG
# ==========================

# Réalise une ORA (Over-Representation Analysis) avec clusterProfiler.
# gene_vec correspond aux gènes différentiellement exprimés à tester.
# universe_vec correspond au background, c'est-à-dire aux gènes testés par
# DESeq2 et possédant une annotation KEGG utilisable.
run_kegg_ora <- function(gene_vec, universe_vec, term2gene_df, label, out_prefix, term2name_df = NULL) {

  gene_vec <- unique(gene_vec)
  universe_vec <- unique(universe_vec)

  message("---- ", label, " ----")
  message("Gènes en entrée : ", length(gene_vec))
  message("Univers         : ", length(universe_vec))
  message("Pathways testés : ", length(unique(term2gene_df[[1]])))

  # Sauvegarde des informations de diagnostic afin de garder une trace
  # du nombre de gènes utilisés pour chaque enrichissement.
  write_csv(
    tibble(
      label = label,
      input_gene_n = length(gene_vec),
      universe_n = length(universe_vec),
      tested_terms_n = length(unique(term2gene_df[[1]]))
    ),
    paste0(out_prefix, "_input_summary.csv")
  )

  if (length(gene_vec) < 3) {
    message("Pas assez de gènes pour enrichissement.")
    return(NULL)
  }

  enr <- tryCatch(
    enricher(
      gene = gene_vec,
      universe = universe_vec,
      TERM2GENE = term2gene_df,
      TERM2NAME = term2name_df,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 5
    ),
    error = function(e) {
      message("Erreur enricher(): ", e$message)
      return(NULL)
    }
  )

  if (is.null(enr)) return(NULL)

  enr_df <- as.data.frame(enr)

  if (nrow(enr_df) == 0) {
    message("Aucun pathway significatif.")
    return(NULL)
  }

  # Filtre complémentaire : conserve les pathways enrichis soutenus par
  # au moins trois gènes, afin de limiter les résultats trop fragiles.
  enr_df <- enr_df %>%
    filter(Count >= 3) %>%
    arrange(p.adjust)

  if (nrow(enr_df) == 0) {
    message("Aucun pathway après filtre Count >= 3.")
    return(NULL)
  }

  write_csv(enr_df, paste0(out_prefix, ".csv"))

  # Dotplot des pathways enrichis.
  p_dot <- dotplot(enr, showCategory = min(20, nrow(enr_df))) +
    ggtitle(label)

  ggsave(
    paste0(out_prefix, "_dotplot.png"),
    plot = p_dot,
    width = 10,
    height = 7,
    dpi = 300
  )

  # Barplot des pathways enrichis.
  p_bar <- barplot(enr, showCategory = min(20, nrow(enr_df))) +
    ggtitle(label)

  ggsave(
    paste0(out_prefix, "_barplot.png"),
    plot = p_bar,
    width = 10,
    height = 7,
    dpi = 300
  )

  message("Pathways significatifs : ", nrow(enr_df))
  return(enr_df)
}

# ==========================
# 6. Analyse par stade expérimental
# ==========================

summary_list <- list()

for (stg in stage_levels) {

  message("===================================")
  message("Stade : ", stg)
  message("===================================")

  de_file <- file.path(deseq_dir, stg, paste0("DE_", stg, "_condition_N_vs_O.csv"))

  if (!file.exists(de_file)) {
    message("Fichier DESeq2 manquant pour ", stg, " : étape ignorée.")
    next
  }

  stage_dir <- file.path(out_dir, stg)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  de <- read_csv(de_file, show_col_types = FALSE)
  colnames(de)[1] <- "gene_id"

  de <- de %>%
    distinct(gene_id, .keep_all = TRUE)

  # Background = gènes réellement testés par DESeq2 au stade considéré.
  # Les gènes sans padj ne sont pas conservés dans l'univers statistique.
  all_genes <- unique(de$gene_id[!is.na(de$padj)])

  # Définition des gènes différentiellement exprimés.
  # Le contraste est N vs O : log2FoldChange positif = plus exprimé en N.
  sig <- de %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= 0.5
    )

  sig_up <- sig %>% filter(log2FoldChange > 0)
  sig_down <- sig %>% filter(log2FoldChange < 0)

  # Restriction du background aux gènes testés par DESeq2 et annotés KEGG.
  ann_kegg_bg <- ann_kegg %>%
    filter(locus_tag %in% all_genes)

  universe_kegg <- unique(ann_kegg_bg$locus_tag)

  # Tableau de diagnostic : utile pour interpréter l'absence éventuelle
  # d'enrichissement significatif à certains stades.
  diag_df <- tibble(
    metric = c(
      "genes_tested_DESeq2",
      "genes_DE_all",
      "genes_DE_up",
      "genes_DE_down",
      "genes_with_KEGG_metabolism_in_background",
      "KEGG_metabolism_pathways_in_background",
      "DE_genes_with_KEGG_metabolism_all",
      "DE_genes_with_KEGG_metabolism_up",
      "DE_genes_with_KEGG_metabolism_down"
    ),
    value = c(
      length(all_genes),
      nrow(sig),
      nrow(sig_up),
      nrow(sig_down),
      length(unique(ann_kegg_bg$locus_tag)),
      length(unique(ann_kegg_bg$KEGG_Pathway)),
      length(intersect(sig$gene_id, universe_kegg)),
      length(intersect(sig_up$gene_id, universe_kegg)),
      length(intersect(sig_down$gene_id, universe_kegg))
    )
  )

  write_csv(diag_df, file.path(stage_dir, paste0("diagnostic_KEGG_metabolism_", stg, ".csv")))

  # Sauvegarde des gènes DE annotés afin de faciliter l'interprétation biologique.
  write_csv(
    sig %>% left_join(ann, by = c("gene_id" = "locus_tag")),
    file.path(stage_dir, paste0("DE_genes_", stg, "_annotated.csv"))
  )

  # Enrichissement des gènes plus exprimés en condition N.
  kegg_up <- run_kegg_ora(
    gene_vec = intersect(sig_up$gene_id, universe_kegg),
    universe_vec = universe_kegg,
    term2gene_df = ann_kegg_bg %>% select(KEGG_Pathway, locus_tag),
    label = paste0("KEGG metabolism - ", stg, " upregulated N vs O"),
    out_prefix = file.path(stage_dir, paste0("KEGG_metabolism_", stg, "_up")),
    term2name_df = kegg_term2name
  )

  # Enrichissement des gènes moins exprimés en condition N.
  kegg_down <- run_kegg_ora(
    gene_vec = intersect(sig_down$gene_id, universe_kegg),
    universe_vec = universe_kegg,
    term2gene_df = ann_kegg_bg %>% select(KEGG_Pathway, locus_tag),
    label = paste0("KEGG metabolism - ", stg, " downregulated N vs O"),
    out_prefix = file.path(stage_dir, paste0("KEGG_metabolism_", stg, "_down")),
    term2name_df = kegg_term2name
  )

  # Résumé du stade courant.
  summary_list[[stg]] <- tibble(
    stage = stg,
    n_genes_tested = length(all_genes),
    n_DE_all = nrow(sig),
    n_DE_up = nrow(sig_up),
    n_DE_down = nrow(sig_down),
    n_KEGG_metabolism_bg_genes = length(universe_kegg),
    n_KEGG_metabolism_bg_pathways = length(unique(ann_kegg_bg$KEGG_Pathway)),
    KEGG_metabolism_up_terms = ifelse(is.null(kegg_up), 0, nrow(kegg_up)),
    KEGG_metabolism_down_terms = ifelse(is.null(kegg_down), 0, nrow(kegg_down))
  )
}

# ==========================
# 7. Résumé final
# ==========================

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  write_csv(summary_df, file.path(out_dir, "KEGG_metabolism_summary_by_stage.csv"))
  print(summary_df, n = nrow(summary_df))
}

message("Analyse terminée.")
message("Résultats écrits dans : ", out_dir)
