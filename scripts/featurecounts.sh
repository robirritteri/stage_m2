#!/bin/bash
#$ -S /bin/bash
#$ -N featurecounts
#$ -q short.q
#$ -V
#$ -cwd
#$ -o /home/rbirritteri/work/logs/sge_featurecounts.out
#$ -e /home/rbirritteri/work/logs/sge_featurecounts.err
#$ -pe thread 4

set -euo pipefail

source /usr/local/genome/miniforge3/etc/profile.d/conda.sh

THREADS="${NSLOTS:-4}"
FEATURECOUNTS="/usr/local/genome/miniforge3/envs/subread-2.0.8/bin/featureCounts"
[[ -x "$FEATURECOUNTS" ]] || { echo "featureCounts introuvable: $FEATURECOUNTS"; exit 1; }

BAM_DIRS=(
  "/home/rbirritteri/work/data/index_mapping/bam"
  "/home/rbirritteri/save/bam"
)
ANN_BASE="/home/rbirritteri/save/annotation"
OUT_DIR="/home/rbirritteri/work/data/counts"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "featureCounts: $FEATURECOUNTS"
echo "THREADS:       $THREADS"
echo "BAM_DIRS:      ${BAM_DIRS[*]}"
echo "ANN_BASE:      $ANN_BASE"
echo "OUT_DIR:       $OUT_DIR"

shopt -s nullglob

# Modifier "bc01..."
for bc in bc01 bc02 bc03 bc04; do
  echo "===== $bc ====="

  GFF="${ANN_BASE}/${bc}_prokka_annot/${bc}.gff"
  [[ -f "$GFF" ]] || { echo "GFF manquant: $GFF"; exit 1; }  

  # Liste complète de BAM pour cette espèce
  ALL_BAMS=()
  for DIR in "${BAM_DIRS[@]}"; do
    ALL_BAMS+=("${DIR}"/*.${bc}.sorted.bam)
  done
  [[ ${#ALL_BAMS[@]} -gt 0 ]] || { echo "Aucun BAM trouvé pour $bc"; continue; }

  # Dédup + tri
  mapfile -t ALL_BAMS_UNIQ < <(printf '%s\n' "${ALL_BAMS[@]}" | sort -u)

  # Sortie
  OUT_PREFIX="${OUT_DIR}/${bc}"
  LOG_FILE="${LOG_DIR}/${bc}.featureCounts.log"

  echo "Nb BAM ($bc): ${#ALL_BAMS_UNIQ[@]}"
  echo "Lancement featureCounts ($bc)..."

  "$FEATURECOUNTS" \
    -T "$THREADS" \
    -F GFF \
    -a "$GFF" \
    -o "${OUT_PREFIX}.counts.txt" \
    -t CDS \
    -g locus_tag \
    -p --countReadPairs -B -C -s 0\
    "${ALL_BAMS_UNIQ[@]}" \
    > "$LOG_FILE" 2>&1

  # Format DESeq2 (garde gène + colonnes de counts)
  grep -v '^#' "${OUT_PREFIX}.counts.txt" | cut -f1,7- > "${OUT_PREFIX}.DESeq2.tsv"

  # Nettoyage noms colonnes: enlève path + .sorted.bam
  awk 'BEGIN{FS="\t"; OFS="\t"}
  NR==1{
    for(i=2;i<=NF;i++){
      gsub(".*/","",$i)
      sub("\\.sorted\\.bam$","",$i)
    }
  }
  {print $0}' "${OUT_PREFIX}.DESeq2.tsv" > "${OUT_PREFIX}.DESeq2.cleaned.tsv"

  echo "OK -> ${OUT_PREFIX}.DESeq2.cleaned.tsv"
  echo "Log -> $LOG_FILE"
  tail -n 12 "$LOG_FILE" || true
done

echo "===== featurecounts fin ====="

# (Optionnel) Génère un samplesheet à partir des colonnes de bc04
# (utile pour DESeq2: condition O/H via impair/pair)
BC_FOR_SAMPLESHEET="bc04"
CLEAN="${OUT_DIR}/${BC_FOR_SAMPLESHEET}.DESeq2.cleaned.tsv"
SAMPLESHEET="${OUT_DIR}/samplesheet.tsv"

if [[ -f "$CLEAN" ]]; then
  echo "Génération samplesheet: $SAMPLESHEET"
  head -n 1 "$CLEAN" \
    | tr '\t' '\n' \
    | tail -n +2 \
    | awk 'BEGIN{OFS="\t"; print "sample","group","num","condition"}
      {
        sample=$1
        # sample type: Prefix-01.bc04  (car on a gardé .bcXX dans le nom)
        gsub(/\..*$/, "", sample)  # enlève .bc04
        split(sample, a, "-")
        group=a[1]
        num=a[2]+0
        if (num <= 4) {
          cond = (num % 2 == 1) ? "O" : "N"
        } else {
          cond = (num % 2 == 1) ? "N" : "O"
        }
        print sample, group, num, cond
      }' > "$SAMPLESHEET"
fi
