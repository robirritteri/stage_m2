#!/bin/bash

# ==========================================================
# Script : parse_eggnog.sh
# Auteur : Romain BIRRITTERI
# Date   : 2026
#
# Description :
#   - Extraction des colonnes fonctionnelles principales du fichier eggNOG-mapper
#   - Génération d'un tableau d'annotation simplifié utilisable dans R
#   - Conservation des informations utiles pour l'interprétation biologique :
#     gène, description, GO, KEGG KO, KEGG Pathway, catégorie COG et enzyme EC
#
# Entrées :
#   - Fichier bc04.emapper.annotations produit par eggNOG-mapper
#
# Sorties :
#   - bc04_annotation.tsv
#
# Usage :
#   bash parse_eggnog.sh
#
# ==========================================================

set -euo pipefail

# ==========================
# Paramètres
# ==========================
# bc04 correspond à Streptococcus thermophilus dans le consortium GazFrom.
BC="bc04"

# ==========================
# Paths
# ==========================
INFILE="/home/rbirritteri/work/data/eggnog/${BC}/${BC}.emapper.annotations"
OUTDIR="/home/rbirritteri/work/data/functional_annotation"
OUTFILE="${OUTDIR}/${BC}_annotation.tsv"

mkdir -p "$OUTDIR"

# ==========================
# Vérification
# ==========================
[[ -f "$INFILE" ]] || { echo "Fichier introuvable: $INFILE"; exit 1; }

# ==========================
# Parsing des annotations eggNOG
# ==========================
# Le fichier eggNOG contient de nombreuses colonnes. Cette étape conserve
# uniquement les colonnes nécessaires aux analyses fonctionnelles réalisées
# ensuite sous R, notamment les annotations KEGG utilisées pour l'enrichissement.
awk -F '\t' 'BEGIN{
    OFS="\t";
    print "locus_tag","gene","description","GO","KEGG_KO","KEGG_Pathway","COG_category","EC"
}
$1 !~ /^#/ {
    locus=$1
    cog=$7
    desc=$8
    gene=$9
    go=$10
    ec=$11
    kegg_ko=$12
    kegg_path=$13
    print locus,gene,desc,go,kegg_ko,kegg_path,cog,ec
}' "$INFILE" > "$OUTFILE"

echo "OK -> $OUTFILE"
