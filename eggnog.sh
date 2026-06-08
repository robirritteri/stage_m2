#!/bin/bash

# Shell à utiliser pour l'exécution du job
#$ -S /bin/bash

# Nom du job
#$ -N eggnog_bc04

# Nom de la queue
#$ -q short.q

# Export de toutes les variables d'environnement
#$ -V

# Lance la commande depuis le répertoire où est lancé le script
#$ -cwd

# Sortie standard
#$ -o /home/rbirritteri/work/logs/sge_eggnog_bc04.out

# Sortie d'erreur
#$ -e /home/rbirritteri/work/logs/sge_eggnog_bc04.err

# Utiliser 8 CPUs
#$ -pe thread 8


# ==========================================================
# Script : eggnog.sh
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Annotation fonctionnelle des protéines de Streptococcus thermophilus (bc04)
#   - Utilisation de eggNOG-mapper à partir du fichier protéique généré par Prokka
#   - Attribution d'orthologues, de descriptions fonctionnelles, de catégories COG,
#     d'identifiants GO, KEGG KO, voies KEGG et enzymes EC
#
# Entrées :
#   - Fichier protéique au format FASTA : bc04.faa
#   - Base de données eggNOG-mapper disponible sur Migale
#
# Sorties :
#   - Fichier d'annotation eggNOG : bc04.emapper.annotations
#   - Fichiers intermédiaires produits par eggNOG-mapper
#
# Usage :
#   qsub eggnog.sh
#
# ==========================================================

set -euo pipefail

# ==========================
# Conda (Migale)
# ==========================
# Initialisation de conda puis activation de l'environnement eggNOG-mapper.
source /usr/local/genome/miniforge3/etc/profile.d/conda.sh
conda activate eggnog-mapper-2.1.13

# Threads : récupère l'allocation SGE.
THREADS="${NSLOTS:-8}"

# ==========================
# Paths
# ==========================
# bc04 correspond à Streptococcus thermophilus dans le consortium GazFrom.
BC="bc04"
FAA="/home/rbirritteri/save/annotation/${BC}_prokka_annot/${BC}.faa"
OUT_DIR="/home/rbirritteri/work/data/eggnog/${BC}"
DATA_DIR="/db/outils/eggnog-mapper-2.1.13"

mkdir -p "$OUT_DIR"

# ==========================
# Vérification des fichiers
# ==========================
# Arrête le script si le fichier protéique ou la base eggNOG sont absents.
[[ -f "$FAA" ]] || { echo "FAA introuvable: $FAA"; exit 1; }
[[ -f "${DATA_DIR}/eggnog.db" ]] || { echo "Base eggNOG introuvable dans: $DATA_DIR"; exit 1; }

# ==========================
# Annotation fonctionnelle eggNOG
# ==========================
# emapper.py recherche des orthologues dans la base eggNOG à partir des séquences
# protéiques, puis associe aux protéines des informations fonctionnelles utiles
# pour l'interprétation biologique et les analyses d'enrichissement.
emapper.py \
  -i "$FAA" \
  --output "$BC" \
  --output_dir "$OUT_DIR" \
  --data_dir "$DATA_DIR" \
  --cpu "$THREADS" \
  --itype proteins \
  --override \
  --target_orthologs all \
  --seed_ortholog_evalue 1e-5

echo "OK -> ${OUT_DIR}/${BC}.emapper.annotations"
