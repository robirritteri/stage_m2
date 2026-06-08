#!/bin/bash

# Shell à utiliser pour l'exécution du job
#$ -S /bin/bash

# Nom du job
#$ -N qc_trimming_qc

# Nom de la queue
#$ -q short.q

# Export de toutes les variables d'environnement
#$ -V

# Sortie standard
#$ -o /home/rbirritteri/work/logs/sge_qc_trim_qc.out

# Sortie d’erreur
#$ -e /home/rbirritteri/work/logs/sge_qc_trim_qc.err

# Lance la commande depuis le répertoire où est lancé le script
#$ -cwd

# Utiliser 4 CPUs
#$ -pe thread 4


# ==========================================================
# Script : qc_trimming_qc.sh
# Auteur : Romain BIRRITTERI
# Date   : 2025-11-19 (modification 2026-02-03)
#
# Description :
#   - QC des données paired-end brutes (FastQC + MultiQC)
#   - Trimming paired-end avec Cutadapt
#   - QC des données trimées (FastQC + MultiQC)
#
# Usage :
#   qsub qc_trimming_qc.sh
#
# ==========================================================

set -euo pipefail

# ==========================
# Conda (Migale)
# ==========================
ENV_FASTQC="fastqc-0.12.1"
ENV_MULTIQC="multiqc-1.27.1"
ENV_CUTADAPT="cutadapt-5.1"

# Initialiser conda (miniforge Migale)
source /usr/local/genome/miniforge3/etc/profile.d/conda.sh

# Threads : récupère l'allocation SGE
THREADS="${NSLOTS:-4}"


# ==========================
# Path
# ==========================
RAW_DIR="/home/rbirritteri/save"
TRIM_DIR="/home/rbirritteri/work/data/trimmed/"
QC_RAW_DIR="/home/rbirritteri/work/data/qc/qc_raw"
QC_TRIM_DIR="/home/rbirritteri/work/data/qc/qc_trim"

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$QC_RAW_DIR" "$QC_TRIM_DIR"


# ==========================
# Paramètres de trimming classiques
# ==========================
QUAL=20
MINLEN=30
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


# ==========================
# 1. QC sur les données brutes
# ==========================
echo "=== [1/4] Contrôle qualité des données brutes ==="
conda run -n "$ENV_FASTQC" fastqc -t "$THREADS" -o "$QC_RAW_DIR" "$RAW_DIR"/*.fastq.gz

echo "=== [1/4] MultiQC des données brutes ==="
conda run -n "$ENV_MULTIQC" multiqc -o "$QC_RAW_DIR" "$QC_RAW_DIR"


# ==========================
# 2. Trimming avec Cutadapt
# ==========================
echo "=== [2/4] Trimming avec Cutadapt ==="

for fq1 in "$RAW_DIR"/*_R1.fastq.gz; do
    fname=$(basename "$fq1")
    base="${fname%_R1.fastq.gz}"
    fq2="$RAW_DIR/${base}_R2.fastq.gz"
    out_fq1="$TRIM_DIR/${base}_R1.trim.fastq.gz"
    out_fq2="$TRIM_DIR/${base}_R2.trim.fastq.gz"

    conda run -n "$ENV_CUTADAPT" cutadapt \
        -j "$THREADS" \
        -q "$QUAL","$QUAL" \
        -m "$MINLEN" \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        -o "$out_fq1" \
        -p "$out_fq2" \
        "$fq1" "$fq2"
done


# ==========================
# 3. QC sur les données trimées
# ==========================
echo "=== [3/4] Contrôle qualité post-trimming ==="
conda run -n "$ENV_FASTQC" fastqc -t "$THREADS" -o "$QC_TRIM_DIR" "$TRIM_DIR"/*.fastq.gz

echo "=== [3/4] MultiQC post-trimming ==="
conda run -n "$ENV_MULTIQC" multiqc -o "$QC_TRIM_DIR" "$QC_TRIM_DIR"


echo "===== qc_trimming_qc fin  ====="
