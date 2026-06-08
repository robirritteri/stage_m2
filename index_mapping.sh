#!/bin/bash

# Shell à utiliser pour l'exécution du job
#$ -S /bin/bash

# Nom du job
#$ -N index_mapping

# Nom de la queue
#$ -q short.q

# Export de toutes les variables d'environnement
#$ -V

# Sortie standard
#$ -o /home/rbirritteri/work/logs/sge_index_mapping.out

# Sortie d'erreur
#$ -e /home/rbirritteri/work/logs/sge_index_mapping.err

# Lance la commande depuis le répertoire où est lancé le script
#$ -cwd

# Utiliser 4 CPUs
#$ -pe thread 4


# ==========================================================
# Script : index_mapping.sh
# Auteur : Romain BIRRITTERI
# Date   : 2026-02-03
#
# Description :
#   - Indexation sur les 4 génomes de références (si nécessaire)
#   - Mapping avec Bowtie2
#   - Passage des SAM en BAM avec samtools
#
# Usage :
#   qsub index_mapping.sh
#
# ==========================================================

set -euo pipefail

# ==========================
# Conda (Migale)
# ==========================
ENV_BOWTIE2="bowtie2-2.5.4"
ENV_SAMTOOLS="samtools-1.21"
BOWTIE2="/usr/local/genome/miniforge3/envs/bowtie2-2.5.4/bin/bowtie2"
BOWTIE2_BUILD="/usr/local/genome/miniforge3/envs/bowtie2-2.5.4/bin/bowtie2-build"
SAMTOOLS="/usr/local/genome/miniforge3/envs/samtools-1.21/bin/samtools"

# Initialiser conda (miniforge Migale)
source /usr/local/genome/miniforge3/etc/profile.d/conda.sh

# Threads : récupère l'allocation SGE
THREADS="${NSLOTS:-4}"


# ==========================
# Paths
# ==========================
SORT_DIR="/home/rbirritteri/work/data/sortmerna"
REF_BASE="/home/rbirritteri/save/annotation"

OUT_DIR="/home/rbirritteri/work/data/index_mapping"
IDX_DIR="${OUT_DIR}/index"
LOG_DIR="${OUT_DIR}/logs"
BAM_DIR="${OUT_DIR}/bam"
SUMMARY="${OUT_DIR}/alignment_summary.tsv"

mkdir -p "$OUT_DIR" "$IDX_DIR" "$LOG_DIR" "$BAM_DIR" 


# ==========================
# Génomes de références
# ==========================
declare -A REF
REF[bc01]="${REF_BASE}/bc01_prokka_annot/bc01.fna"
REF[bc02]="${REF_BASE}/bc02_prokka_annot/bc02.fna"
REF[bc03]="${REF_BASE}/bc03_prokka_annot/bc03.fna"
REF[bc04]="${REF_BASE}/bc04_prokka_annot/bc04.fna"


# ==========================
# 1. Indexation Bowtie2
# ==========================
echo "Construction des index Bowtie2 (si nécessaire)"
for bc in bc01 bc02 bc03 bc04; do
  if [[ -f "${IDX_DIR}/${bc}.1.bt2" || -f "${IDX_DIR}/${bc}.1.bt2l" ]]; then
    echo "Index déjà présent pour $bc"
  else
    echo "Construction index $bc"
    "$BOWTIE2_BUILD" --threads "$THREADS" "${REF[$bc]}" "${IDX_DIR}/${bc}"
  fi
done


# ==========================
# 2. Alignement
# ==========================
echo -e "sample\tgenome\toverall_alignment_rate" > "$SUMMARY"

R1_FILES=("${SORT_DIR}"/Fp3*.non_rRNA.fq.gz)

for r1 in "${R1_FILES[@]}"; do
  sample=$(basename "$r1" .non_rRNA.fq.gz)

  echo "Sample: $sample"

  for bc in bc01 bc02 bc03 bc04; do
    echo "Alignement sur $bc"

    bam_out="${BAM_DIR}/${sample}.${bc}.sorted.bam"
    log_out="${LOG_DIR}/${sample}.${bc}.bowtie2.log"

    if ! "$BOWTIE2" \
      -x "${IDX_DIR}/${bc}" \
      --interleaved  "$r1" \
      -p "$THREADS" \
      --very-sensitive \
      --phred33 \
      -S /dev/stdout \
      2> "$log_out" \
     | "$SAMTOOLS" view -@ "$THREADS" -b - \
     | "$SAMTOOLS" sort -@ "$THREADS" -o "$bam_out" - ; then
      echo "Échec alignement sur $bc" >&2
      echo -e "${sample}\t${bc}\tFAIL" >> "$SUMMARY"
      continue
    fi

    "$SAMTOOLS" index "$bam_out"

    rate=$(grep -E "overall alignment rate" "$log_out" | awk '{print $1}' | tail -n 1 || true)
    [[ -n "$rate" ]] || rate="NA"
    echo -e "${sample}\t${bc}\t${rate}" >> "$SUMMARY"
  done
done

echo "BAM:     $BAM_DIR"
echo "Logs:    $LOG_DIR"
echo "Résumé:  $SUMMARY"

echo "===== Index_mapping fin  ====="
