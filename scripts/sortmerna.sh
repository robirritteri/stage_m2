#!/bin/bash

# Shell à utiliser pour l'exécution du job
#$ -S /bin/bash

# Nom du job
#$ -N sortmerna

# Nom de la queue
#$ -q short.q

# Export de toutes les variables d'environnement
#$ -V

# Sortie standard
#$ -o /home/rbirritteri/work/logs/sge_sortmerna.out

# Sortie d’erreur
#$ -e /home/rbirritteri/work/logs/sge_sortmerna.err

# Lance la commande depuis le répertoire où est lancé le script
#$ -cwd

# Utiliser 8 CPUs
#$ -pe thread 8


# ==========================================================
# Script : sortmerna.sh
# Auteur : Romain BIRRITTERI
# Date   : 2026-02-06
#
# Description :
#   - Déplétion des reads rRNA sur données RNA-seq paired-end
#   - Filtrage basé sur les bases SILVA bactériennes (16S / 23S) et RFAM bactériennes (5S / 5.8S)
#   - Conservation des reads non-rRNA pour l’analyse transcriptomique
#
# Usage :
#   qsub sortmerna.sh
#
# ==========================================================

set -euo pipefail

# ==========================
# Conda (Migale)
# ==========================
SORTMERNA="/usr/local/genome/miniforge3/envs/sortmerna-4.3.7/bin/sortmerna"
ENV_FASTQC="fastqc-0.12.1"
ENV_MULTIQC="multiqc-1.27.1"

# Initialiser conda (miniforge Migale)
source /usr/local/genome/miniforge3/etc/profile.d/conda.sh

# Threads : récupère l'allocation SGE
THREADS="${NSLOTS:-8}"

# ==========================
# Paths
# ==========================
IN_DIR="/home/rbirritteri/work/data/trimmed/"
OUT_DIR="/home/rbirritteri/work/data/sortmerna/"
LOG_DIR="${OUT_DIR}/logs"
WORKDIR_BASE="${OUT_DIR}/workdir_db"
QC_SORT_DIR="/home/rbirritteri/work/data/qc/qc_sort"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$WORKDIR_BASE" "$QC_SORT_DIR"  

# ==========================
# Bases rRNA Migale
# ==========================
DB_DIR="/db/outils/sortmerna"

REF1="${DB_DIR}/silva-bac-16s-id90.fasta"
REF2="${DB_DIR}/silva-bac-23s-id98.fasta"
REF3="${DB_DIR}/rfam-5s-database-id98.fasta"
REF4="${DB_DIR}/rfam-5.8s-database-id98.fasta"

# ==========================
# Vérification
# ==========================
[[ -x "$SORTMERNA" ]] || { echo "sortmerna introuvable: $SORTMERNA"; exit 1; }

for f in "$REF1" "$REF2" "$REF3" "$REF4"; do
  [[ -f "$f" ]] || { echo "Fichier DB/index manquant: $f"; exit 1; }
done

shopt -s nullglob
R1_FILES=("${IN_DIR}"/*_R1.trim.fastq.gz)
[[ ${#R1_FILES[@]} -gt 0 ]] || { echo "Aucun *_R1.trim.fastq.gz dans $IN_DIR"; exit 1; }

# ==========================
# Filtrage rRNA
# ==========================
 for r1 in "${R1_FILES[@]}"; do
  sample=$(basename "$r1" _R1.trim.fastq.gz)
  r2="${IN_DIR}/${sample}_R2.trim.fastq.gz"
  
  [[ -f "$r2" ]] || { echo "R2 manquant pour $sample (attendu: $r2)"; exit 1; }

  echo "==> Sample: $sample"
  # Préfixes de sortie
  aligned_prefix="${OUT_DIR}/${sample}.rRNA"
  other_prefix="${OUT_DIR}/${sample}.non_rRNA"
  
  WORKDIR="${WORKDIR_BASE}/${sample}"
  mkdir -p "$WORKDIR"
  
  "$SORTMERNA" \
    --ref "${REF1}" \
    --ref "${REF2}" \
    --ref "${REF3}" \
    --ref "${REF4}" \
    --reads "$r1" \
    --reads "$r2" \
    --paired_in \
    --fastx \
    --other "$other_prefix" \
    --aligned "$aligned_prefix" \
    --workdir "$WORKDIR" \
    --threads "$THREADS" \
    > "${LOG_DIR}/${sample}.sortmerna.log" 2>&1

done

echo "Terminé. Résultats dans: $OUT_DIR"
echo "non-rRNA : ${OUT_DIR}/*.non_rRNA*"
echo "rRNA     : ${OUT_DIR}/*.rRNA*"


echo "=== [1/2] Contrôle qualité post-déplétion ==="
conda run -n "$ENV_FASTQC" fastqc -t "$THREADS" -o "$QC_SORT_DIR" "$OUT_DIR"/*.non_rRNA.*

echo "=== [2/2] MultiQC post-déplétion ==="
conda run -n "$ENV_MULTIQC" multiqc -o "$QC_SORT_DIR" "$QC_SORT_DIR"


echo "===== Contrôle qualité post-déplétion fin  ====="
